/*
 * cfe_serialize.c — CFE state checkpoint/restore implementation
 *
 * Implements the ngen serialization protocol (see ngen_utilities.h).
 * Engine call sequence:
 *   Save:    SetValue(create) -> GetValue(size) -> GetValue(state) -> SetValue(free)
 *   Restore: SetValue(state, payload)
 *
 * See cfe_serialize.h for the developer checklist when changing serialized state.
 *
 * Layout version 1 byte buffer:
 *
 * Array variables are always preceded by their element count so the
 * deserializer can validate before reading.
 *
 *   [4]   uint32_t  layout_version = 1
 *   [4]   int       flags (bit 0 = DSBM active)
 *   [8]   double    soil_storage_m
 *   [4]   int       n_disc D
 *  [8*D]  double[D] soil_discrete_storage_theta
 *   [8]   double    gw_storage_m
 *   [4]   int       n_nash_subsurface S
 *  [8*S]  double[S] nash_subsurface_storage_m
 *   [4]   int       n_giuh G
 *  [8*G]  double[G] giuh_queue_m
 *   [8]   double    skin_temperature_k
 *   [8]   double    upper_soil_temperature_k
 *   [8]   double    estimated_annual_air_temp_k
 *   [8]   double    air_temp_time_integral_k_s
 *   [8]   double    accumulated_time_s
 *   [4]   int       pet_initialized
 *  Total: 80 + 8*(D+S+G) bytes  (default: 80 + 8*(4+2+G))
 */

#include <stdio.h>
#include <stdlib.h>
#include "cfe_serialize.h"
#include "ngen_utilities.h"
#include "bmi.h"

static int serialization_flags(const CFE_Model_Context *ctx) {
    int flags = 0;
    if (ctx->options.simulate_discrete_soil_moisture) flags |= CFE_SERIALIZATION_FLAG_DSBM;
    return flags;
}

size_t cfe_serialize_total_bytes(const CFE_Model_Context *ctx) {
    return CFE_SERIALIZATION_FIXED_BYTES
        + (size_t)ctx->parameters.giuh_num_ordinates * sizeof(double);
}

void cfe_serialize_create(CFE_Model_Context *ctx) {
    if (ctx->serialized_state != NULL)
        free(ctx->serialized_state);

    ctx->serialized_size = cfe_serialize_total_bytes(ctx);
    ctx->serialized_state = (char *)malloc(ctx->serialized_size);

    char *p = ctx->serialized_state;

    p = ser_write_u32(p, CFE_SERIALIZATION_LAYOUT_VERSION);
    p = ser_write_i32(p, serialization_flags(ctx));

    p = ser_write_f64(p, ctx->state.soil_storage_m);

    p = ser_write_i32(p, NDISC);
    p = ser_write_f64_array(p, ctx->state.soil_discrete_storage_theta, NDISC);

    p = ser_write_f64(p, ctx->state.gw_storage_m);

    p = ser_write_i32(p, MAX_NUM_SUBSURFACE_NASH_CASCADE);
    p = ser_write_f64_array(p, ctx->state.nash_subsurface_storage_m,
                            MAX_NUM_SUBSURFACE_NASH_CASCADE);

    p = ser_write_i32(p, ctx->parameters.giuh_num_ordinates);
    p = ser_write_f64_array(p, ctx->state.giuh_queue_m,
                        ctx->parameters.giuh_num_ordinates);

    /* PET temperature state */
    p = ser_write_f64(p, ctx->state.pet_temperature_state.skin_temperature_k);
    p = ser_write_f64(p, ctx->state.pet_temperature_state.upper_soil_temperature_k);
    p = ser_write_f64(p, ctx->state.pet_temperature_state.estimated_annual_air_temperature_k);
    p = ser_write_f64(p, ctx->state.pet_temperature_state.air_temperature_time_integral_k_s);
    p = ser_write_f64(p, ctx->state.pet_temperature_state.accumulated_time_s);
    ser_write_i32(p, ctx->state.pet_temperature_state.initialized);
}

int cfe_serialize_deserialize(CFE_Model_Context *ctx, const char *src) {
    const char *p = src;

    uint32_t version;
    p = ser_read_u32(p, &version);
    if (version != CFE_SERIALIZATION_LAYOUT_VERSION) {
        fprintf(stderr, "cfe_serialize_deserialize: checkpoint layout version %u "
                "does not match this build (expected %d). "
                "The checkpoint was created by an incompatible CFE version.\n",
                version, CFE_SERIALIZATION_LAYOUT_VERSION);
        return BMI_FAILURE;
    }

    int flags;
    p = ser_read_i32(p, &flags);
    if (flags != serialization_flags(ctx)) {
        fprintf(stderr, "cfe_serialize_deserialize: checkpoint flags 0x%x "
                "do not match current config flags 0x%x. "
                "Check that DSBM settings match between save and restore configs.\n",
                flags, serialization_flags(ctx));
        return BMI_FAILURE;
    }

    p = ser_read_f64(p, &ctx->state.soil_storage_m);

    int n_disc;
    p = ser_read_i32(p, &n_disc);
    if (n_disc != NDISC) {
        fprintf(stderr, "cfe_serialize_deserialize: checkpoint has n_disc=%d "
                "but this build was compiled with NDISC=%d. "
                "The checkpoint was created by a build with a different soil discretization.\n",
                n_disc, NDISC);
        return BMI_FAILURE;
    }
    p = ser_read_f64_array(p, ctx->state.soil_discrete_storage_theta, n_disc);

    p = ser_read_f64(p, &ctx->state.gw_storage_m);

    int n_nash;
    p = ser_read_i32(p, &n_nash);
    if (n_nash != MAX_NUM_SUBSURFACE_NASH_CASCADE) {
        fprintf(stderr, "cfe_serialize_deserialize: checkpoint has n_nash_subsurface=%d "
                "but this build was compiled with MAX_NUM_SUBSURFACE_NASH_CASCADE=%d. "
                "The checkpoint was created by a build with a different Nash cascade size.\n",
                n_nash, MAX_NUM_SUBSURFACE_NASH_CASCADE);
        return BMI_FAILURE;
    }
    p = ser_read_f64_array(p, ctx->state.nash_subsurface_storage_m, n_nash);

    int n_giuh;
    p = ser_read_i32(p, &n_giuh);
    if (n_giuh != ctx->parameters.giuh_num_ordinates) {
        fprintf(stderr, "cfe_serialize_deserialize: checkpoint has n_giuh=%d "
                "but the current config has giuh_num_ordinates=%d. "
                "The checkpoint was created with a different GIUH configuration.\n",
                n_giuh, ctx->parameters.giuh_num_ordinates);
        return BMI_FAILURE;
    }
    p = ser_read_f64_array(p, ctx->state.giuh_queue_m, n_giuh);

    /* PET temperature state */
    p = ser_read_f64(p, &ctx->state.pet_temperature_state.skin_temperature_k);
    p = ser_read_f64(p, &ctx->state.pet_temperature_state.upper_soil_temperature_k);
    p = ser_read_f64(p, &ctx->state.pet_temperature_state.estimated_annual_air_temperature_k);
    p = ser_read_f64(p, &ctx->state.pet_temperature_state.air_temperature_time_integral_k_s);
    p = ser_read_f64(p, &ctx->state.pet_temperature_state.accumulated_time_s);
    int pet_init;
    p = ser_read_i32(p, &pet_init);
    ctx->state.pet_temperature_state.initialized = pet_init;

    /* Sync DSBM derived state from restored theta values.
     * cfe() reads soil_state_in.theta_in before DSBM updates it, so stale
     * values after a restore would corrupt rainfall partitioning. */
    if (ctx->options.simulate_discrete_soil_moisture) {
        ctx->state.soil_state_in.total_storage_m = 0.0;
        for (int i = 0; i < NDISC; i++) {
            ctx->state.soil_state_in.theta_in[i] = ctx->state.soil_discrete_storage_theta[i];
            ctx->state.soil_state_in.total_storage_m +=
                ctx->state.soil_state_in.theta_in[i] * ctx->state.soil_geometry.dz_m[i];
        }
    }

    return BMI_SUCCESS;
}

void cfe_serialize_free(CFE_Model_Context *ctx) {
    if (ctx->serialized_state != NULL) {
        free(ctx->serialized_state);
        ctx->serialized_state = NULL;
    }
    ctx->serialized_size = 0;
}
