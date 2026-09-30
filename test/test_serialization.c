/*
 * test_serialization.c — unit tests for cfe_serialize.c
 *
 * Tests the serialization layer directly (not through BMI dispatch).
 * Uses BMI only to initialize and advance the model to build up state.
 *
 * Usage: test_serialization <test_name>
 *   test_serialize_create       — buffer allocation, header contents, free
 *   test_deserialize_failures   — rejects corrupted buffers
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "bmi.h"
#include "bmi_cfe.h"
#include "cfe_serialize.h"
#include "ngen_utilities.h"

#define PASS 0
#define FAIL 1

static Bmi bmi_storage;

static int init_and_advance(Bmi *m, const char *cfg, int steps) {
    register_bmi_cfe(m);
    if (m->initialize(m, cfg) != BMI_SUCCESS) return FAIL;
    double rain = 0.005, pet = 0.0;
    for (int t = 0; t < steps; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        if (m->update(m) != BMI_SUCCESS) return FAIL;
    }
    return PASS;
}

static int test_serialize_create(const char *cfg) {
    Bmi *m = &bmi_storage;
    if (init_and_advance(m, cfg, 3) != PASS) {
        printf("  FAIL: init_and_advance\n");
        return FAIL;
    }

    CFE_Model_Context *ctx = (CFE_Model_Context *)m->data;

    /* --- Create --- */
    cfe_serialize_create(ctx);

    size_t expected_size = cfe_serialize_total_bytes(ctx);
    if (ctx->serialized_state == NULL) {
        printf("  FAIL: serialized_state is NULL after create\n");
        return FAIL;
    }
    if (ctx->serialized_size != expected_size) {
        printf("  FAIL: serialized_size=%zu, expected=%zu\n",
               ctx->serialized_size, expected_size);
        return FAIL;
    }

    /* Walk the buffer and verify header + array counts */
    const char *p = ctx->serialized_state;

    uint32_t version;
    p = ser_read_u32(p, &version);
    if (version != CFE_SERIALIZATION_LAYOUT_VERSION) {
        printf("  FAIL: buffer version=%u, expected=%d\n",
               version, CFE_SERIALIZATION_LAYOUT_VERSION);
        return FAIL;
    }

    int flags;
    p = ser_read_i32(p, &flags);
    int expected_flags = ctx->options.simulate_discrete_soil_moisture
                         ? CFE_SERIALIZATION_FLAG_DSBM : 0;
    if (flags != expected_flags) {
        printf("  FAIL: buffer flags=0x%x, expected=0x%x\n", flags, expected_flags);
        return FAIL;
    }

    /* soil_storage_m scalar */
    double soil_val;
    p = ser_read_f64(p, &soil_val);
    if (soil_val != ctx->state.soil_storage_m) {
        printf("  FAIL: buffer soil_storage_m=%.15e, expected=%.15e\n",
               soil_val, ctx->state.soil_storage_m);
        return FAIL;
    }

    /* n_disc count */
    int n_disc;
    p = ser_read_i32(p, &n_disc);
    if (n_disc != NDISC) {
        printf("  FAIL: buffer n_disc=%d, expected NDISC=%d\n", n_disc, NDISC);
        return FAIL;
    }

    /* skip theta array */
    double dummy_theta[NDISC];
    p = ser_read_f64_array(p, dummy_theta, NDISC);

    /* gw_storage_m scalar */
    double gw_val;
    p = ser_read_f64(p, &gw_val);
    if (gw_val != ctx->state.gw_storage_m) {
        printf("  FAIL: buffer gw_storage_m=%.15e, expected=%.15e\n",
               gw_val, ctx->state.gw_storage_m);
        return FAIL;
    }

    /* n_nash count */
    int n_nash;
    p = ser_read_i32(p, &n_nash);
    if (n_nash != MAX_NUM_SUBSURFACE_NASH_CASCADE) {
        printf("  FAIL: buffer n_nash=%d, expected=%d\n",
               n_nash, MAX_NUM_SUBSURFACE_NASH_CASCADE);
        return FAIL;
    }

    /* skip nash array */
    double dummy_nash[MAX_NUM_SUBSURFACE_NASH_CASCADE];
    p = ser_read_f64_array(p, dummy_nash, MAX_NUM_SUBSURFACE_NASH_CASCADE);

    /* n_giuh count */
    int n_giuh;
    p = ser_read_i32(p, &n_giuh);
    if (n_giuh != ctx->parameters.giuh_num_ordinates) {
        printf("  FAIL: buffer n_giuh=%d, expected=%d\n",
               n_giuh, ctx->parameters.giuh_num_ordinates);
        return FAIL;
    }

    /* --- Double-create (should free + reallocate without leak) --- */
    cfe_serialize_create(ctx);
    if (ctx->serialized_state == NULL || ctx->serialized_size != expected_size) {
        printf("  FAIL: double-create\n");
        return FAIL;
    }

    /* --- Free --- */
    cfe_serialize_free(ctx);
    if (ctx->serialized_state != NULL) {
        printf("  FAIL: serialized_state not NULL after free\n");
        return FAIL;
    }
    if (ctx->serialized_size != 0) {
        printf("  FAIL: serialized_size=%zu after free, expected 0\n",
               ctx->serialized_size);
        return FAIL;
    }

    /* Double-free should be safe */
    cfe_serialize_free(ctx);

    printf("  PASS: buffer %zu bytes, header validated, "
           "double-create and double-free safe\n", expected_size);
    m->finalize(m);
    return PASS;
}

static int test_deserialize_failures(const char *cfg) {
    Bmi *m = &bmi_storage;
    if (init_and_advance(m, cfg, 3) != PASS) {
        printf("  FAIL: init_and_advance\n");
        return FAIL;
    }

    CFE_Model_Context *ctx = (CFE_Model_Context *)m->data;

    /* Create a valid reference buffer */
    cfe_serialize_create(ctx);
    size_t buf_size = ctx->serialized_size;
    char *good_buf = (char *)malloc(buf_size);
    memcpy(good_buf, ctx->serialized_state, buf_size);
    cfe_serialize_free(ctx);

    char *bad_buf = (char *)malloc(buf_size);
    int result = PASS;

    /* Compute field offsets by walking the layout */
    size_t off_version = 0;
    size_t off_flags   = off_version + sizeof(uint32_t);
    size_t off_soil    = off_flags + sizeof(int);
    size_t off_n_disc  = off_soil + sizeof(double);
    size_t off_theta   = off_n_disc + sizeof(int);
    size_t off_gw      = off_theta + NDISC * sizeof(double);
    size_t off_n_nash  = off_gw + sizeof(double);
    size_t off_nash    = off_n_nash + sizeof(int);
    size_t off_n_giuh  = off_nash + MAX_NUM_SUBSURFACE_NASH_CASCADE * sizeof(double);

    /* Case 1: wrong layout version */
    memcpy(bad_buf, good_buf, buf_size);
    uint32_t bad_version = 99;
    memcpy(bad_buf + off_version, &bad_version, sizeof(bad_version));
    if (cfe_serialize_deserialize(ctx, bad_buf) != BMI_FAILURE) {
        printf("  FAIL: expected rejection for wrong layout version\n");
        result = FAIL;
    }

    /* Case 2: wrong flags */
    memcpy(bad_buf, good_buf, buf_size);
    int bad_flags = 0xFF;
    memcpy(bad_buf + off_flags, &bad_flags, sizeof(bad_flags));
    if (cfe_serialize_deserialize(ctx, bad_buf) != BMI_FAILURE) {
        printf("  FAIL: expected rejection for wrong flags\n");
        result = FAIL;
    }

    /* Case 3: wrong n_disc */
    memcpy(bad_buf, good_buf, buf_size);
    int bad_n_disc = NDISC + 1;
    memcpy(bad_buf + off_n_disc, &bad_n_disc, sizeof(bad_n_disc));
    if (cfe_serialize_deserialize(ctx, bad_buf) != BMI_FAILURE) {
        printf("  FAIL: expected rejection for wrong n_disc\n");
        result = FAIL;
    }

    /* Case 4: wrong n_nash */
    memcpy(bad_buf, good_buf, buf_size);
    int bad_n_nash = MAX_NUM_SUBSURFACE_NASH_CASCADE + 1;
    memcpy(bad_buf + off_n_nash, &bad_n_nash, sizeof(bad_n_nash));
    if (cfe_serialize_deserialize(ctx, bad_buf) != BMI_FAILURE) {
        printf("  FAIL: expected rejection for wrong n_nash\n");
        result = FAIL;
    }

    /* Case 5: wrong n_giuh */
    memcpy(bad_buf, good_buf, buf_size);
    int bad_n_giuh = ctx->parameters.giuh_num_ordinates + 1;
    memcpy(bad_buf + off_n_giuh, &bad_n_giuh, sizeof(bad_n_giuh));
    if (cfe_serialize_deserialize(ctx, bad_buf) != BMI_FAILURE) {
        printf("  FAIL: expected rejection for wrong n_giuh\n");
        result = FAIL;
    }

    /* Sanity: valid buffer should succeed */
    if (cfe_serialize_deserialize(ctx, good_buf) != BMI_SUCCESS) {
        printf("  FAIL: valid buffer unexpectedly rejected\n");
        result = FAIL;
    }

    if (result == PASS)
        printf("  PASS: all 5 corruption cases rejected, valid buffer accepted\n");

    free(good_buf);
    free(bad_buf);
    m->finalize(m);
    return result;
}

int main(int argc, const char *argv[]) {
    if (argc < 2) {
        printf("Usage: test_serialization <test_name> [config]\n");
        printf("  test_serialize_create\n");
        printf("  test_deserialize_failures\n");
        return FAIL;
    }

    const char *cfg = (argc >= 3) ? argv[2] : TEST_SERIALIZATION_CONFIG;

    if (strcmp(argv[1], "test_serialize_create") == 0)
        return test_serialize_create(cfg);
    if (strcmp(argv[1], "test_deserialize_failures") == 0)
        return test_deserialize_failures(cfg);

    printf("Unknown test: %s\n", argv[1]);
    return FAIL;
}
