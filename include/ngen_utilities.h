/*
Author: Nels Frazier
Copyright (C) 2025 Lynker
------------------------------------------------------------------------
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
------------------------------------------------------------------------

Version 0.1
Protocol utilities for use with the ngen model engine
*/
#include <stdint.h>
#include <stddef.h>

#define NGEN_MASS_IN "ngen::mass_in"
#define NGEN_MASS_OUT "ngen::mass_out"
#define NGEN_MASS_STORED "ngen::mass_stored"
#define NGEN_MASS_LEAKED "ngen::mass_leaked"

#define MASS_BALANCE_VAR_NAME_COUNT 4

#include <bmi.h>
#include <string.h>

//These are the bare minimum required to interface with ngen
static const char *mass_balance_var_names[MASS_BALANCE_VAR_NAME_COUNT] = { NGEN_MASS_IN, NGEN_MASS_OUT, NGEN_MASS_STORED, NGEN_MASS_LEAKED };
static const char *mass_balance_var_types[MASS_BALANCE_VAR_NAME_COUNT] = { "double", "double", "double", "double" };
static const int   mass_balance_var_item_count[MASS_BALANCE_VAR_NAME_COUNT] = { 1, 1, 1, 1};
static const char *mass_balance_var_units[MASS_BALANCE_VAR_NAME_COUNT] = { "m", "m", "m", "m" };
//These may be useful in the future, but aren't strictly necessary for the protocol at the moment
// static const char *mass_balance_var_grids[MASS_BALANCE_VAR_NAME_COUNT] = { 0, 0, 0 };
// static const char *mass_balance_var_locations[MASS_BALANCE_VAR_NAME_COUNT] = { "node", "node", "node" };

static inline int get_mass_balance_var_type(const char *name, char *type)
{
    for (int i = 0; i < MASS_BALANCE_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, mass_balance_var_names[i]) == 0) {
            strncpy(type, mass_balance_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    return BMI_FAILURE;
}

static inline int get_mass_balance_item_count(const char *name)
{
    for (int i = 0; i < MASS_BALANCE_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, mass_balance_var_names[i]) == 0) {
            return mass_balance_var_item_count[i];
        }
    }
    return -1; // Default to -1 if not found
}

static inline int get_mass_balance_unit(const char *name, char* unit)
{
    for (int i = 0; i < MASS_BALANCE_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, mass_balance_var_names[i]) == 0) {
            strncpy(unit, mass_balance_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    return BMI_FAILURE;
}

/*
 * Serialization Protocol — opt-in BMI state checkpoint/restore
 *
 * These reserved variables allow a BMI model to participate in the ngen
 * state serialization workflow (see NOAA-OWP/ngen#957).  They are NOT
 * advertised via GetInputVarNames/GetOutputVarNames but must be
 * resolvable by name through GetVarType, GetVarUnits, etc.
 *
 * Engine save sequence:
 *   SetValue(create) -> GetValue(size) -> GetValue(state, buf) -> write -> SetValue(free)
 * Engine restore:
 *   SetValue(state, payload)
 */

#define NGEN_SERIALIZATION_CREATE "ngen::serialization_create"
#define NGEN_SERIALIZATION_FREE   "ngen::serialization_free"
#define NGEN_SERIALIZATION_SIZE   "ngen::serialization_size"
#define NGEN_SERIALIZATION_STATE  "ngen::serialization_state"

#define SERIALIZATION_VAR_NAME_COUNT 4

static const char *serialization_var_names[SERIALIZATION_VAR_NAME_COUNT] = {
    NGEN_SERIALIZATION_CREATE, NGEN_SERIALIZATION_FREE,
    NGEN_SERIALIZATION_SIZE, NGEN_SERIALIZATION_STATE
};
static const char *serialization_var_types[SERIALIZATION_VAR_NAME_COUNT] = {
    "int", "int", "int", "char"
};
static const int serialization_var_item_count[SERIALIZATION_VAR_NAME_COUNT] = {
    1, 1, 1, 0
};
static const char *serialization_var_units[SERIALIZATION_VAR_NAME_COUNT] = {
    "ngen::trigger", "ngen::trigger", "bytes", "ngen::opaque"
};

static inline int get_serialization_var_type(const char *name, char *type)
{
    for (int i = 0; i < SERIALIZATION_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, serialization_var_names[i]) == 0) {
            strncpy(type, serialization_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    return BMI_FAILURE;
}

static inline int get_serialization_item_count(const char *name)
{
    for (int i = 0; i < SERIALIZATION_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, serialization_var_names[i]) == 0) {
            return serialization_var_item_count[i];
        }
    }
    return -1;
}

static inline int get_serialization_unit(const char *name, char *unit)
{
    for (int i = 0; i < SERIALIZATION_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, serialization_var_names[i]) == 0) {
            strncpy(unit, serialization_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    return BMI_FAILURE;
}

static inline int is_serialization_var(const char *name)
{
    for (int i = 0; i < SERIALIZATION_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, serialization_var_names[i]) == 0) {
            return 1;
        }
    }
    return 0;
}

/* Get itemsize for serialization variables.
 * Triggers (create/free) have no stored value -> BMI_FAILURE.
 * size -> sizeof(int), state -> sizeof(char). */
static inline int get_serialization_itemsize(const char *name, int *size)
{
    if (strcmp(name, NGEN_SERIALIZATION_CREATE) == 0 ||
        strcmp(name, NGEN_SERIALIZATION_FREE)   == 0) {
        *size = 0;
        return BMI_FAILURE;
    }
    if (strcmp(name, NGEN_SERIALIZATION_SIZE) == 0) {
        *size = sizeof(int);
        return BMI_SUCCESS;
    }
    if (strcmp(name, NGEN_SERIALIZATION_STATE) == 0) {
        *size = sizeof(char);
        return BMI_SUCCESS;
    }
    return BMI_FAILURE;
}

/* Get nbytes for serialization variables.
 * Triggers -> BMI_FAILURE (no stored value).
 * size -> sizeof(int).
 * state -> requires current buffer size from model, returns BMI_FAILURE
 *          here (caller must handle state separately with model context). */
static inline int get_serialization_nbytes(const char *name, int *nbytes)
{
    if (strcmp(name, NGEN_SERIALIZATION_CREATE) == 0 ||
        strcmp(name, NGEN_SERIALIZATION_FREE)   == 0) {
        *nbytes = 0;
        return BMI_FAILURE;
    }
    if (strcmp(name, NGEN_SERIALIZATION_SIZE) == 0) {
        *nbytes = sizeof(int);
        return BMI_SUCCESS;
    }
    /* NGEN_SERIALIZATION_STATE: caller must handle — size is dynamic */
    return BMI_FAILURE;
}

/* ================================================================== */
/*  Buffer-cursor helpers for serialization protocol                   */
/* ================================================================== */
/* Model-agnostic primitives: memcpy + pointer advance in one call.
 * Any BMI module implementing the serialization protocol can use these. */

static inline char *ser_write_u32(char *dest, uint32_t val) {
    memcpy(dest, &val, sizeof(val));
    return dest + sizeof(val);
}

static inline char *ser_write_i32(char *dest, int val) {
    memcpy(dest, &val, sizeof(val));
    return dest + sizeof(val);
}

static inline char *ser_write_f64(char *dest, double val) {
    memcpy(dest, &val, sizeof(val));
    return dest + sizeof(val);
}

static inline char *ser_write_f64_array(char *dest, const double *src, int n) {
    size_t bytes = (size_t)n * sizeof(double);
    memcpy(dest, src, bytes);
    return dest + bytes;
}

static inline const char *ser_read_u32(const char *src, uint32_t *val) {
    memcpy(val, src, sizeof(*val));
    return src + sizeof(*val);
}

static inline const char *ser_read_i32(const char *src, int *val) {
    memcpy(val, src, sizeof(*val));
    return src + sizeof(*val);
}

static inline const char *ser_read_f64(const char *src, double *val) {
    memcpy(val, src, sizeof(*val));
    return src + sizeof(*val);
}

static inline const char *ser_read_f64_array(const char *src, double *dest, int n) {
    size_t bytes = (size_t)n * sizeof(double);
    memcpy(dest, src, bytes);
    return src + bytes;
}