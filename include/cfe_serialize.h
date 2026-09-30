#ifndef CFE_SERIALIZE_H
#define CFE_SERIALIZE_H

/*
 * CFE state serialization — ngen checkpoint/restore protocol
 *
 * Only evolving state variables are serialized — fields marked [SERIALIZED]
 * in cfe_state_struct (cfe_types.h).  Parameters, options, and lookup tables
 * are reconstructed from the config file during Initialize.
 *
 * Convention: array variables are always preceded by their element count
 * so the deserializer can validate before reading.
 *
 * All serialized arrays in CFE are fixed-size, declared with compile-time
 * constants (NDISC, MAX_NUM_SUBSURFACE_NASH_CASCADE) or config values
 * (giuh_num_ordinates).  This means the array data is inline in the
 * struct — no pointers to chase — and a flat memcpy is sufficient.
 * If a future field uses a dynamically allocated array (e.g. double*),
 * the serializer must dereference the pointer and write the pointed-to
 * data, and the deserializer must allocate before reading into it.
 * The count-then-array convention still applies, but the buffer size
 * can no longer be computed at compile time from the macros below.
 *
 * ADDING OR REMOVING STATE — developer checklist:
 *   1. Mark the field [SERIALIZED] in cfe_state_struct (cfe_types.h)
 *   2. Add a row to the byte layout diagram in cfe_serialize.c
 *   3. Add cursor-helper calls in cfe_serialize_create() and
 *      cfe_serialize_deserialize() (cfe_serialize.c)
 *      For arrays: write the count (ser_write_i32) then the array
 *   4. Update CFE_SERIALIZATION_FIXED_BYTES below (or the variable-length
 *      section in cfe_serialize_total_bytes())
 *   5. Bump CFE_SERIALIZATION_LAYOUT_VERSION so old buffers are rejected
 *   6. If the new field has derived state that the model reads before the
 *      first timestep updates it, sync that derived state at the end of
 *      cfe_serialize_deserialize() (see the soil_state_in.theta_in sync)
 *   7. Add the field to test_serialization_round_trip and
 *      test_serialization_round_trip_dsbm in test_bmi_model.c
 */

#include <stdint.h>
#include <stddef.h>
#include "cfe_types.h"

#define CFE_SERIALIZATION_LAYOUT_VERSION 1
#define CFE_SERIALIZATION_FLAG_DSBM      0x1
#define CFE_SERIALIZATION_HEADER_BYTES  (sizeof(uint32_t) + sizeof(int))
#define CFE_SERIALIZATION_FIXED_BYTES   (CFE_SERIALIZATION_HEADER_BYTES \
    + sizeof(double)                     /* soil_storage_m */ \
    + sizeof(int)                        /* n_disc count */ \
    + NDISC * sizeof(double)             /* soil_discrete_storage_theta */ \
    + sizeof(double)                     /* gw_storage_m */ \
    + sizeof(int)                        /* n_nash_subsurface count */ \
    + MAX_NUM_SUBSURFACE_NASH_CASCADE * sizeof(double)  /* nash_subsurface */ \
    + sizeof(int)                        /* n_giuh count */ \
    + 5 * sizeof(double)                 /* PET temperature state (5 doubles) */ \
    + sizeof(int))                       /* pet_initialized */

void   cfe_serialize_create(CFE_Model_Context *ctx);
int    cfe_serialize_deserialize(CFE_Model_Context *ctx, const char *src);
void   cfe_serialize_free(CFE_Model_Context *ctx);
size_t cfe_serialize_total_bytes(const CFE_Model_Context *ctx);

#endif
