/*
 * bmi_cfe.c — CFE v3 BMI implementation
 *
 * Based on cfe3-project/src/cfe_bmi.c by FLO, adapted to the CSDMS BMI
 * function structure.  Model state lives in CFE_Model_Context stored in
 * self->data (accessed via the CONTEXT macro).
 *
 *  - Initialize delegates to cfe_context_create_from_config()
 *  - Update delegates to cfe_context_update()
 *  - Finalize delegates to cfe_context_destroy()
 *  - get_value_ptr supports all variables + ngen mass balance protocol
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include "bmi.h"
#include "bmi_cfe.h"
#include "cfe_context.h"
#include "cfe_helpers.h"
#include "parser_helpers.h"
#include "ngen_utilities.h"
#include "cfe_serialize.h"

/* ------------------------------------------------------------------ */
/* Cast helper — extract CFE_Model_Context from the BMI data pointer  */
/* ------------------------------------------------------------------ */
#define CONTEXT(self) ((CFE_Model_Context*)(self)->data)

/* ================================================================== */
/*  Variable name tables                                               */
/* ================================================================== */

/* --- inputs (model forcing only) ---
 * verbosity is accessible via set_value/get_value/get_value_ptr but is
 * not advertised as a BMI input variable.
 * The CSDMS BMI standard has no "string" type, and ngen rejects non-numeric inputs. */
static const char* input_var_names[] = {
    "rainfall_depth_m",
    "et_potential_m",
    "ice_fraction",
    "day_of_year",
    "DLWRF_surface",
    "DSWRF_surface",
    "PRES_surface",
    "SPFH_2maboveground",
    "TMP_2maboveground",
    "UGRD_10maboveground",
    "VGRD_10maboveground",
    "param_catchment_vegetated_fraction",
    "bare_soil_rsurf_exp"
};
static const int INPUT_VAR_NAME_COUNT = 13;

/* --- outputs --- */
static const char* output_var_names[] = {
    /* Primary fluxes (per timestep) */
    "discharge_m",                              /*  0 */
    "surface_runoff_m",                         /*  1 */
    "lateral_flow_m",                           /*  2 */
    "baseflow_m",                               /*  3 */
    "actual_et_m",                              /*  4 */
    "vol_balance_residual_m",                   /*  5 */

    /* State scalars (checkpointing) */
    "state_soil_storage_m",                     /*  6 */
    "state_gw_storage_m",                       /*  7 */
    "state_current_timestep",                   /*  8 */

    /* State arrays (checkpointing) */
    "state_soil_moisture_theta",                /*  9  NDISC elements */
    "state_nash_subsurface_storage",            /* 10  2 elements */
    "state_giuh_queue",                         /* 11  num_giuh elements */

    /* Config / parameters for interpretation */
    "config_simulate_discrete_soil_moisture",   /* 12 */
    "param_catchment_area_km2",                 /* 13 */
    "param_soil_depth_m",                       /* 14 */
    "param_soil_porosity",                      /* 15 */

    /* Per-timestep volume balance */
    "timestep_storage_start_m",                 /* 16 */
    "timestep_input_m",                         /* 17 */
    "timestep_output_m",                        /* 18 */
    "timestep_storage_end_m",                   /* 19 */

    /* Additional flux outputs */
    "potential_et_m",                           /* 20 */
    "giuh_outflow_m",                           /* 21 */
    "soil_to_gw_percolation_flux_m",            /* 22 */

    /* Per-layer DSBM soil moisture (scalar aliases for ngen CSV output) */
    "soil_moisture_theta_1",                    /* 23 */
    "soil_moisture_theta_2",                    /* 24 */
    "soil_moisture_theta_3",                    /* 25 */
    "soil_moisture_theta_4",                    /* 26 */

    /* Detailed flux breakdown */
    "bare_soil_evaporation_m",                  /* 27 */
    "impervious_runoff_m",                      /* 28 */
    "pervious_runoff_m",                        /* 29 */
    "surface_routed_to_outlet_m",               /* 30 */
    "lateral_flow_generated_m"                  /* 31 */
};
static const int OUTPUT_VAR_NAME_COUNT = 32;

/* --- calibration parameters (get_value / set_value / get_value_ptr) --- */
static const char* param_var_names[] = {
    "soil_effective_porosity",                   /*  0  dimensionless */
    "soil_saturated_hydraulic_conductivity",     /*  1  cm/h */
    "soil_percolation_rate_limiter",             /*  2  0-1 */
    "soil_Clapp_Hornberger_b",                  /*  3  dimensionless */
    "soil_lateral_flow_K",                       /*  4  per h */
    "subsurface_nash_K",                         /*  5  per h */
    "gw_discharge_coefficient",                  /*  6  m/s */
    "gw_discharge_exponent",                     /*  7  dimensionless */
    "gw_max_storage_m",                          /*  8  m */
    "soil_saturated_capillary_head",             /*  9  cm */
    "soil_field_capacity_fraction",              /* 10  Pcap/Patm */
    "Xinanjiang_inflection_a",                   /* 11  0-1 */
    "Xinanjiang_shape_b",                        /* 12  dimensionless */
    "Xinanjiang_shape_x",                        /* 13  dimensionless */
    "Priestley_Taylor_alpha",                    /* 14  dimensionless */
    "soil_ice_imperv_threshold"                  /* 15  dimensionless */
};
static const int PARAM_VAR_NAME_COUNT = 16;

/* Resolve a v3 parameter name to a pointer into ctx->parameters.
 * Returns NULL if unrecognized. */
static double* param_field_ptr(CFE_Model_Context *ctx, const char *name) {
    cfe_parameters_struct *p = &ctx->parameters;

    if (strcmp(name, "soil_effective_porosity") == 0)               return &p->effective_porosity;
    if (strcmp(name, "soil_saturated_hydraulic_conductivity") == 0) return &p->ksat_m_per_s;
    if (strcmp(name, "soil_percolation_rate_limiter") == 0)         return &p->soil_to_gw_percolation_rate_limiter_0_1;
    if (strcmp(name, "soil_Clapp_Hornberger_b") == 0)              return &p->soil_b;
    if (strcmp(name, "soil_lateral_flow_K") == 0)                   return &p->soil_k_lateral_per_h;
    if (strcmp(name, "subsurface_nash_K") == 0)                     return &p->nash_subsurface_K_per_h;
    if (strcmp(name, "gw_discharge_coefficient") == 0)              return &p->gw_discharge_coeff_m_per_timestep;
    if (strcmp(name, "gw_discharge_exponent") == 0)                 return &p->gw_discharge_exponent;
    if (strcmp(name, "gw_max_storage_m") == 0)                      return &p->gw_max_storage_m;
    if (strcmp(name, "soil_saturated_capillary_head") == 0)         return &p->sat_capillary_head_m;
    if (strcmp(name, "soil_field_capacity_fraction") == 0)          return &p->field_capacity_Pcap_over_Patm;
    if (strcmp(name, "Xinanjiang_inflection_a") == 0)               return &p->xj_tension_inflection_point;
    if (strcmp(name, "Xinanjiang_shape_b") == 0)                    return &p->xj_tension_b;
    if (strcmp(name, "Xinanjiang_shape_x") == 0)                    return &p->xj_free_b;
    if (strcmp(name, "Priestley_Taylor_alpha") == 0)                return &p->alpha_pt;
    if (strcmp(name, "soil_ice_imperv_threshold") == 0)             return &p->soil_ice_imperv_threshold;

    return NULL;
}

/* ================================================================== */
/*  Lifecycle: Initialize / Update / Finalize                          */
/* ================================================================== */

static int Initialize(Bmi *self, const char *cfg_file) {
    if (cfg_file == NULL || self == NULL) return BMI_FAILURE;

    CFE_Model_Context* ctx = NULL;
    int result = cfe_context_create_from_config(cfg_file, &ctx);
    if (result != 0 || ctx == NULL) return BMI_FAILURE;

    self->data = (void*)ctx;

    /* initial ngen mass balance: cumulative_vol starts at initial total storage
     * so that mass_in always represents total mass in the domain */
    double init_storage = calculate_total_storage(ctx);
    ctx->volbal.cumulative_vol   = init_storage;
    ctx->volbal.volume_in_domain = init_storage;
    ctx->volbal.leakage          = 0.0;

    return BMI_SUCCESS;
}

static int Update(Bmi *self) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;

    CFE_Model_Context *ctx = CONTEXT(self);

    /* BMI inputs arrive as rates (m/s); convert to depth (m) for the timestep */
    double dt = (double)ctx->options.time_step_seconds;
    ctx->forcing.rainfall_depth_m *= dt;
    ctx->forcing.et_potential_m   *= dt;

    /* accumulate input BEFORE the step so mass_in is consistent on failure */
    ctx->volbal.cumulative_vol += ctx->forcing.rainfall_depth_m;

    int result = cfe_context_update(ctx);

    /* update remaining protocol fields regardless of success/failure
     * so the caller can inspect mass balance state on BMI_FAILURE */
    ctx->volbal.volume_in_domain = ctx->timestep_storage_end_m;
    ctx->vol_balance_residual_m = ctx->volbal.volstart + ctx->volbal.volin
                                - ctx->volbal.volout   - ctx->volbal.volend;

    if (result != 0) return BMI_FAILURE;
    return BMI_SUCCESS;
}

static int Update_until(Bmi *self, double then) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;

    double current_time = 0.0;
    int dt_seconds;
    cfe_context_get_time_step_seconds(CONTEXT(self), &dt_seconds);
    int step;
    cfe_context_get_current_step(CONTEXT(self), &step);
    current_time = (double)step * (double)dt_seconds;

    while (current_time < then) {
        int result = Update(self);
        if (result != BMI_SUCCESS) return result;
        cfe_context_get_current_step(CONTEXT(self), &step);
        current_time = (double)step * (double)dt_seconds;
    }
    return BMI_SUCCESS;
}

static int Finalize(Bmi *self) {
    if (CONTEXT(self) != NULL) {
        cfe_context_destroy(CONTEXT(self));
        self->data = NULL;
    }
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Component / exchange item metadata                                 */
/* ================================================================== */

static int Get_component_name(Bmi *self, char *name) {
    strcpy(name, "CFE - Conceptual Functional Equivalent");
    return BMI_SUCCESS;
}

static int needs_aorc_forcing(Bmi *self) {
    if (CONTEXT(self) == NULL) return 0;
    return CONTEXT(self)->options.enable_ET_Priestley_Taylor ||
           CONTEXT(self)->options.simulate_soil_evaporation;
}

static int has_internal_day_of_year(Bmi *self) {
    return CONTEXT(self) != NULL && CONTEXT(self)->options.epoch_start_seconds > 0.0;
}

static int is_aorc_only_input(const char *name) {
    return strcmp(name, "day_of_year") == 0 ||
           strcmp(name, "DLWRF_surface") == 0 ||
           strcmp(name, "DSWRF_surface") == 0 ||
           strcmp(name, "PRES_surface") == 0 ||
           strcmp(name, "SPFH_2maboveground") == 0 ||
           strcmp(name, "TMP_2maboveground") == 0 ||
           strcmp(name, "UGRD_10maboveground") == 0 ||
           strcmp(name, "VGRD_10maboveground") == 0;
}

static int should_skip_input(Bmi *self, const char *name) {
    if (!needs_aorc_forcing(self) && is_aorc_only_input(name))
        return 1;
    if (has_internal_day_of_year(self) && strcmp(name, "day_of_year") == 0)
        return 1;
    return 0;
}

static int count_skipped_inputs(Bmi *self) {
    int skip = 0;
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++)
        if (should_skip_input(self, input_var_names[i])) skip++;
    return skip;
}

static int Get_input_item_count(Bmi *self, int *count) {
    *count = INPUT_VAR_NAME_COUNT - count_skipped_inputs(self);
    return BMI_SUCCESS;
}

static int Get_output_item_count(Bmi *self, int *count) {
    *count = OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}

static int Get_input_var_names(Bmi *self, char **names) {
    int j = 0;
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++) {
        if (should_skip_input(self, input_var_names[i]))
            continue;
        strcpy(names[j++], input_var_names[i]);
    }
    return BMI_SUCCESS;
}

static int Get_output_var_names(Bmi *self, char **names) {
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++)
        strcpy(names[i], output_var_names[i]);
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Variable metadata                                                  */
/* ================================================================== */

static int Get_var_grid(Bmi *self, const char *name, int *grid) {
    if      (strcmp(name, "state_soil_moisture_theta")    == 0) *grid = 1;
    else if (strcmp(name, "state_nash_subsurface_storage")== 0) *grid = 2;
    else if (strcmp(name, "state_giuh_queue")             == 0) *grid = 3;
    else *grid = 0;
    return BMI_SUCCESS;
}

static int Get_var_type(Bmi *self, const char *name, char *type) {
    /* Serialization protocol variables */
    if (get_serialization_var_type(name, type) == BMI_SUCCESS)
        return BMI_SUCCESS;

    if (strcmp(name, "verbosity") == 0 ||
        strcmp(name, "state_current_timestep") == 0 ||
        strcmp(name, "config_simulate_discrete_soil_moisture") == 0 ||
        strcmp(name, "day_of_year") == 0) {
        strcpy(type, "int");
        return BMI_SUCCESS;
    }
    /* Check output and input variable names (all double) */
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++)
        if (strcmp(name, output_var_names[i]) == 0) { strcpy(type, "double"); return BMI_SUCCESS; }
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++)
        if (strcmp(name, input_var_names[i]) == 0) { strcpy(type, "double"); return BMI_SUCCESS; }
    /* Check calibration parameter names (all double) */
    if (CONTEXT(self) && param_field_ptr(CONTEXT(self), name) != NULL) {
        strcpy(type, "double");
        return BMI_SUCCESS;
    }
    return BMI_FAILURE;
}

static int Get_var_units(Bmi *self, const char *name, char *units) {
    /* Serialization protocol variables */
    if (get_serialization_unit(name, units) == BMI_SUCCESS)
        return BMI_SUCCESS;

    if (strcmp(name, "rainfall_depth_m")         == 0 ||
        strcmp(name, "et_potential_m")           == 0) {
        strcpy(units, "m s-1");
        return BMI_SUCCESS;
    }
    if (strcmp(name, "discharge_m")              == 0 ||
        strcmp(name, "surface_runoff_m")         == 0 ||
        strcmp(name, "lateral_flow_m")           == 0 ||
        strcmp(name, "baseflow_m")               == 0 ||
        strcmp(name, "actual_et_m")              == 0 ||
        strcmp(name, "vol_balance_residual_m")   == 0 ||
        strcmp(name, "state_soil_storage_m")     == 0 ||
        strcmp(name, "state_gw_storage_m")       == 0 ||
        strcmp(name, "param_soil_depth_m")       == 0 ||
        strcmp(name, "timestep_storage_start_m") == 0 ||
        strcmp(name, "timestep_input_m")         == 0 ||
        strcmp(name, "timestep_output_m")        == 0 ||
        strcmp(name, "timestep_storage_end_m")   == 0 ||
        strcmp(name, "potential_et_m")          == 0 ||
        strcmp(name, "giuh_outflow_m")          == 0 ||
        strcmp(name, "soil_to_gw_percolation_flux_m") == 0 ||
        strcmp(name, "bare_soil_evaporation_m") == 0 ||
        strcmp(name, "impervious_runoff_m") == 0 ||
        strcmp(name, "pervious_runoff_m") == 0 ||
        strcmp(name, "surface_routed_to_outlet_m") == 0 ||
        strcmp(name, "lateral_flow_generated_m") == 0) {
        strcpy(units, "m");
    }
    else if (strcmp(name, "state_soil_moisture_theta") == 0 ||
             strcmp(name, "param_soil_porosity") == 0 ||
             strcmp(name, "soil_moisture_theta_1") == 0 ||
             strcmp(name, "soil_moisture_theta_2") == 0 ||
             strcmp(name, "soil_moisture_theta_3") == 0 ||
             strcmp(name, "soil_moisture_theta_4") == 0) {
        strcpy(units, "-");
    }
    else if (strcmp(name, "state_nash_subsurface_storage") == 0 ||
             strcmp(name, "state_giuh_queue")              == 0) {
        strcpy(units, "m");
    }
    else if (strcmp(name, "param_catchment_area_km2") == 0) {
        strcpy(units, "km2");
    }
    else if (strcmp(name, "verbosity") == 0 ||
             strcmp(name, "state_current_timestep") == 0 ||
             strcmp(name, "config_simulate_discrete_soil_moisture") == 0 ||
             strcmp(name, "day_of_year") == 0) {
        strcpy(units, "1");
    }
    /* --- AORC forcing variables --- */
    else if (strcmp(name, "DLWRF_surface") == 0 ||
             strcmp(name, "DSWRF_surface") == 0) {
        strcpy(units, "W m-2");
    }
    else if (strcmp(name, "PRES_surface") == 0) {
        strcpy(units, "Pa");
    }
    else if (strcmp(name, "SPFH_2maboveground") == 0) {
        strcpy(units, "kg kg-1");
    }
    else if (strcmp(name, "TMP_2maboveground") == 0) {
        strcpy(units, "K");
    }
    else if (strcmp(name, "UGRD_10maboveground") == 0 ||
             strcmp(name, "VGRD_10maboveground") == 0) {
        strcpy(units, "m s-1");
    }
    /* --- calibration parameter units (internal representation) --- */
    else if (strcmp(name, "soil_effective_porosity") == 0 ||
             strcmp(name, "soil_field_capacity_fraction") == 0 ||
             strcmp(name, "soil_percolation_rate_limiter") == 0 ||
             strcmp(name, "Xinanjiang_inflection_a") == 0 ||
             strcmp(name, "soil_ice_imperv_threshold") == 0 ||
             strcmp(name, "ice_fraction") == 0) {
        strcpy(units, "-");  /* dimensionless fractions (V/V or 0-1) */
    }
    else if (strcmp(name, "soil_saturated_hydraulic_conductivity") == 0) {
        strcpy(units, "cm h-1");
    }
    else if (strcmp(name, "soil_saturated_capillary_head") == 0) {
        strcpy(units, "cm");
    }
    else if (strcmp(name, "gw_discharge_coefficient") == 0) {
        strcpy(units, "m s-1");
    }
    else if (strcmp(name, "soil_lateral_flow_K") == 0 ||
             strcmp(name, "subsurface_nash_K") == 0) {
        strcpy(units, "h-1");
    }
    else if (strcmp(name, "gw_max_storage_m") == 0) {
        strcpy(units, "m");
    }
    else if (strcmp(name, "soil_Clapp_Hornberger_b") == 0 ||
             strcmp(name, "gw_discharge_exponent") == 0 ||
             strcmp(name, "Xinanjiang_shape_b") == 0 ||
             strcmp(name, "Xinanjiang_shape_x") == 0 ||
             strcmp(name, "Priestley_Taylor_alpha") == 0 ||
             strcmp(name, "bare_soil_rsurf_exp") == 0) {
        strcpy(units, "-");  /* dimensionless exponents and coefficients */
    }
    else if (strcmp(name, "param_catchment_vegetated_fraction") == 0) {
        strcpy(units, "-");
    }
    else {
        return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

static int Get_var_itemsize(Bmi *self, const char *name, int *size) {
    /* Serialization protocol variables */
    if (is_serialization_var(name))
        return get_serialization_itemsize(name, size);

    if (strcmp(name, "verbosity") == 0 ||
        strcmp(name, "state_current_timestep") == 0 ||
        strcmp(name, "config_simulate_discrete_soil_moisture") == 0 ||
        strcmp(name, "day_of_year") == 0) {
        *size = sizeof(int);
    }
    else {
        /* all remaining recognized variables (outputs, inputs, params) are double */
        char type[BMI_MAX_TYPE_NAME];
        if (Get_var_type(self, name, type) != BMI_SUCCESS) return BMI_FAILURE;
        if (strcmp(type, "double") == 0) *size = sizeof(double);
        else if (strcmp(type, "int") == 0) *size = sizeof(int);
        else return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

static int Get_var_nbytes(Bmi *self, const char *name, int *nbytes) {
    /* Serialization protocol variables */
    if (is_serialization_var(name)) {
        if (get_serialization_nbytes(name, nbytes) == BMI_SUCCESS)
            return BMI_SUCCESS;
        /* state: dynamic size from model context */
        if (strcmp(name, NGEN_SERIALIZATION_STATE) == 0) {
            *nbytes = CONTEXT(self) ? (int)CONTEXT(self)->serialized_size : 0;
            return BMI_SUCCESS;
        }
        return BMI_FAILURE;  /* triggers */
    }

    if (strcmp(name, "verbosity") == 0 ||
        strcmp(name, "state_current_timestep") == 0 ||
        strcmp(name, "config_simulate_discrete_soil_moisture") == 0 ||
        strcmp(name, "day_of_year") == 0) {
        *nbytes = sizeof(int);
    }
    else if (strcmp(name, "state_soil_moisture_theta") == 0) {
        *nbytes = NDISC * sizeof(double);
    }
    else if (strcmp(name, "state_nash_subsurface_storage") == 0) {
        *nbytes = 2 * sizeof(double);
    }
    else if (strcmp(name, "state_giuh_queue") == 0) {
        if (CONTEXT(self) && CONTEXT(self)->options.surface_routing_scheme == SURF_ROUTE_GIUH)
            *nbytes = CONTEXT(self)->parameters.giuh_num_ordinates * sizeof(double);
        else
            *nbytes = 0;
    }
    else {
        /* remaining recognized variables (scalar outputs, inputs, params) are single doubles */
        int itemsize;
        if (Get_var_itemsize(self, name, &itemsize) != BMI_SUCCESS) return BMI_FAILURE;
        *nbytes = itemsize;
    }
    return BMI_SUCCESS;
}

static int Get_var_location(Bmi *self, const char *name, char *location) {
    /* Serialization protocol vars have no spatial semantics */
    if (is_serialization_var(name)) return BMI_FAILURE;

    /* Verify this is a recognized variable before returning location */
    char type[BMI_MAX_TYPE_NAME];
    if (Get_var_type(self, name, type) != BMI_SUCCESS) return BMI_FAILURE;
    strcpy(location, "node");
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Time information                                                   */
/* ================================================================== */

static int Get_current_time(Bmi *self, double *time) {
    if (CONTEXT(self) == NULL) { *time = 0.0; return BMI_SUCCESS; }
    int step;
    cfe_context_get_current_step(CONTEXT(self), &step);
    int dt_s;
    cfe_context_get_time_step_seconds(CONTEXT(self), &dt_s);
    *time = (double)step * (double)dt_s;
    return BMI_SUCCESS;
}

static int Get_start_time(Bmi *self, double *time) {
    *time = 0.0;
    return BMI_SUCCESS;
}

static int Get_end_time(Bmi *self, double *time) {
    if (CONTEXT(self) != NULL && CONTEXT(self)->options.num_timesteps > 0) {
        int dt_s;
        cfe_context_get_time_step_seconds(CONTEXT(self), &dt_s);
        *time = (double)CONTEXT(self)->options.num_timesteps * (double)dt_s;
    } else {
        *time = (double)FLT_MAX;  /* unknown — forcings arrive via BMI */
    }
    return BMI_SUCCESS;
}

static int Get_time_units(Bmi *self, char *units) {
    strcpy(units, "s");
    return BMI_SUCCESS;
}

static int Get_time_step(Bmi *self, double *dt) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;
    int dt_s;
    cfe_context_get_time_step_seconds(CONTEXT(self), &dt_s);
    *dt = (double)dt_s;
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Get value / Get value ptr / Get value at indices                   */
/* ================================================================== */

/* forward declarations — Get_value delegates to Get_value_ptr */
static int Get_value_ptr(Bmi *self, const char *name, void **dest);
static int Get_var_nbytes(Bmi *self, const char *name, int *nbytes);

static int Get_value(Bmi *self, const char *name, void *dest) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;

    /* Serialization protocol */
    if (strcmp(name, NGEN_SERIALIZATION_SIZE) == 0) {
        /* BMI GetVarNbytes uses int; serialized_size is size_t.
         * Safe as long as serialized buffers stay under INT_MAX (~2 GB). */
        *(int*)dest = (int)CONTEXT(self)->serialized_size;
        return BMI_SUCCESS;
    }
    if (strcmp(name, NGEN_SERIALIZATION_STATE) == 0) {
        if (CONTEXT(self)->serialized_state && CONTEXT(self)->serialized_size > 0)
            memcpy(dest, CONTEXT(self)->serialized_state, CONTEXT(self)->serialized_size);
        return BMI_SUCCESS;
    }

    /* Arrays need element-by-element copy from the source array */
    if (strcmp(name, "state_soil_moisture_theta") == 0) {
        double *d = (double*)dest;
        for (int i = 0; i < NDISC; i++)
            d[i] = CONTEXT(self)->state.soil_discrete_storage_theta[i];
        return BMI_SUCCESS;
    }
    if (strcmp(name, "state_nash_subsurface_storage") == 0) {
        double *d = (double*)dest;
        for (int i = 0; i < 2; i++)
            d[i] = CONTEXT(self)->state.nash_subsurface_storage_m[i];
        return BMI_SUCCESS;
    }
    if (strcmp(name, "state_giuh_queue") == 0) {
        double *d = (double*)dest;
        for (int i = 0; i < CONTEXT(self)->parameters.giuh_num_ordinates; i++)
            d[i] = CONTEXT(self)->state.giuh_queue_m[i];
        return BMI_SUCCESS;
    }
    /* Params with BMI-boundary unit conversion (internal SI → user-facing) */
    if (strcmp(name, "soil_saturated_hydraulic_conductivity") == 0) {
        *(double*)dest = m_per_s_to_cm_per_h(CONTEXT(self)->parameters.ksat_m_per_s);
        return BMI_SUCCESS;
    }
    if (strcmp(name, "soil_saturated_capillary_head") == 0) {
        *(double*)dest = m_to_cm(CONTEXT(self)->parameters.sat_capillary_head_m);
        return BMI_SUCCESS;
    }

    /* All other variables: delegate through get_value_ptr */
    void *ptr = NULL;
    if (Get_value_ptr(self, name, &ptr) != BMI_SUCCESS || ptr == NULL)
        return BMI_FAILURE;

    /* Determine size to copy */
    int nbytes = 0;
    if (Get_var_nbytes(self, name, &nbytes) != BMI_SUCCESS)
        return BMI_FAILURE;

    memcpy(dest, ptr, nbytes);
    return BMI_SUCCESS;
}

static int Get_value_ptr(Bmi *self, const char *name, void **dest) {
    CFE_Model_Context *ctx = CONTEXT(self);
    if (ctx == NULL) return BMI_FAILURE;

    /* --- ngen serialization protocol --- */
    if (strcmp(name, NGEN_SERIALIZATION_SIZE) == 0)  { *dest = &ctx->serialized_size;  return BMI_SUCCESS; }
    if (strcmp(name, NGEN_SERIALIZATION_STATE) == 0) { *dest = ctx->serialized_state;  return BMI_SUCCESS; }
    /* Triggers (create/free) have no stored value — fall through to BMI_FAILURE */

    /* --- ngen mass balance protocol --- */
    if (strcmp(name, NGEN_MASS_IN) == 0)     { *dest = &ctx->volbal.cumulative_vol;   return BMI_SUCCESS; }
    if (strcmp(name, NGEN_MASS_OUT) == 0)    { *dest = &ctx->volbal.volout;           return BMI_SUCCESS; }
    if (strcmp(name, NGEN_MASS_STORED) == 0) { *dest = &ctx->volbal.volume_in_domain; return BMI_SUCCESS; }
    if (strcmp(name, NGEN_MASS_LEAKED) == 0) { *dest = &ctx->volbal.leakage;          return BMI_SUCCESS; }

    /* --- output scalars --- */
    if (strcmp(name, "discharge_m") == 0)              { *dest = &ctx->last_outputs.qout_m;                    return BMI_SUCCESS; }
    if (strcmp(name, "surface_runoff_m") == 0)         { *dest = &ctx->last_outputs.surface_runoff_generated_m; return BMI_SUCCESS; }
    if (strcmp(name, "lateral_flow_m") == 0)           { *dest = &ctx->last_outputs.lateral_flow_m;            return BMI_SUCCESS; }
    if (strcmp(name, "baseflow_m") == 0)               { *dest = &ctx->last_outputs.baseflow_m;                return BMI_SUCCESS; }
    if (strcmp(name, "actual_et_m") == 0)              { *dest = &ctx->last_outputs.actual_et_m;               return BMI_SUCCESS; }
    if (strcmp(name, "state_soil_storage_m") == 0)     { *dest = &ctx->state.soil_storage_m;                   return BMI_SUCCESS; }
    if (strcmp(name, "state_gw_storage_m") == 0)       { *dest = &ctx->state.gw_storage_m;                     return BMI_SUCCESS; }
    if (strcmp(name, "state_current_timestep") == 0)   { *dest = &ctx->state.current_time_step;                return BMI_SUCCESS; }
    if (strcmp(name, "config_simulate_discrete_soil_moisture") == 0) { *dest = &ctx->options.simulate_discrete_soil_moisture; return BMI_SUCCESS; }
    if (strcmp(name, "param_catchment_area_km2") == 0) { *dest = &ctx->parameters.catchment_area_km2;         return BMI_SUCCESS; }
    if (strcmp(name, "param_soil_depth_m") == 0)       { *dest = &ctx->parameters.soil_depth_m;                return BMI_SUCCESS; }
    if (strcmp(name, "param_soil_porosity") == 0)      { *dest = &ctx->parameters.effective_porosity;          return BMI_SUCCESS; }
    if (strcmp(name, "timestep_storage_start_m") == 0) { *dest = &ctx->timestep_storage_start_m;               return BMI_SUCCESS; }
    if (strcmp(name, "timestep_input_m") == 0)         { *dest = &ctx->timestep_input_m;                       return BMI_SUCCESS; }
    if (strcmp(name, "timestep_output_m") == 0)        { *dest = &ctx->timestep_output_m;                      return BMI_SUCCESS; }
    if (strcmp(name, "timestep_storage_end_m") == 0)   { *dest = &ctx->timestep_storage_end_m;                 return BMI_SUCCESS; }
    if (strcmp(name, "vol_balance_residual_m") == 0) { *dest = &ctx->vol_balance_residual_m;                return BMI_SUCCESS; }
    if (strcmp(name, "potential_et_m") == 0)           { *dest = &ctx->last_outputs.potential_et_m;             return BMI_SUCCESS; }
    if (strcmp(name, "giuh_outflow_m") == 0)           { *dest = &ctx->last_outputs.giuh_outflow_m;            return BMI_SUCCESS; }
    if (strcmp(name, "soil_to_gw_percolation_flux_m") == 0) { *dest = &ctx->last_outputs.soil_to_gw_percolation_flux_m; return BMI_SUCCESS; }
    if (strcmp(name, "bare_soil_evaporation_m") == 0)      { *dest = &ctx->last_outputs.bare_soil_evaporation_m;      return BMI_SUCCESS; }
    if (strcmp(name, "impervious_runoff_m") == 0)           { *dest = &ctx->last_outputs.impervious_runoff_m;          return BMI_SUCCESS; }
    if (strcmp(name, "pervious_runoff_m") == 0)             { *dest = &ctx->last_outputs.pervious_runoff_m;            return BMI_SUCCESS; }
    if (strcmp(name, "surface_routed_to_outlet_m") == 0)    { *dest = &ctx->last_outputs.surface_routed_to_outlet_m;   return BMI_SUCCESS; }
    if (strcmp(name, "lateral_flow_generated_m") == 0)      { *dest = &ctx->last_outputs.lateral_flow_generated_m;     return BMI_SUCCESS; }

    /* --- per-layer DSBM soil moisture scalars --- */
    if (strcmp(name, "soil_moisture_theta_1") == 0) { *dest = &ctx->state.soil_discrete_storage_theta[0]; return BMI_SUCCESS; }
    if (strcmp(name, "soil_moisture_theta_2") == 0) { *dest = &ctx->state.soil_discrete_storage_theta[1]; return BMI_SUCCESS; }
    if (strcmp(name, "soil_moisture_theta_3") == 0) { *dest = &ctx->state.soil_discrete_storage_theta[2]; return BMI_SUCCESS; }
    if (strcmp(name, "soil_moisture_theta_4") == 0) { *dest = &ctx->state.soil_discrete_storage_theta[3]; return BMI_SUCCESS; }

    /* --- state arrays (BMI output variables for checkpointing/hotstart) --- */
    if (strcmp(name, "state_soil_moisture_theta") == 0)     { *dest = ctx->state.soil_discrete_storage_theta;  return BMI_SUCCESS; }
    if (strcmp(name, "state_nash_subsurface_storage") == 0) { *dest = ctx->state.nash_subsurface_storage_m;    return BMI_SUCCESS; }
    if (strcmp(name, "state_giuh_queue") == 0)              { *dest = ctx->state.giuh_queue_m;                 return BMI_SUCCESS; }

    /* --- inputs --- */
    if (strcmp(name, "rainfall_depth_m") == 0)  { *dest = &ctx->forcing.rainfall_depth_m;          return BMI_SUCCESS; }
    if (strcmp(name, "et_potential_m") == 0)    { *dest = &ctx->forcing.et_potential_m;             return BMI_SUCCESS; }
    if (strcmp(name, "ice_fraction") == 0)      { *dest = &ctx->forcing.ice_fraction;              return BMI_SUCCESS; }
    if (strcmp(name, "day_of_year") == 0)       { *dest = &ctx->forcing.day_of_year;               return BMI_SUCCESS; }
    if (strcmp(name, "DLWRF_surface") == 0)     { *dest = &ctx->forcing.DLWRF_surface;             return BMI_SUCCESS; }
    if (strcmp(name, "DSWRF_surface") == 0)     { *dest = &ctx->forcing.DSWRF_surface;             return BMI_SUCCESS; }
    if (strcmp(name, "PRES_surface") == 0)      { *dest = &ctx->forcing.PRES_surface;              return BMI_SUCCESS; }
    if (strcmp(name, "SPFH_2maboveground") == 0){ *dest = &ctx->forcing.SPFH_2maboveground;       return BMI_SUCCESS; }
    if (strcmp(name, "TMP_2maboveground") == 0) { *dest = &ctx->forcing.TMP_2maboveground;         return BMI_SUCCESS; }
    if (strcmp(name, "UGRD_10maboveground") == 0){*dest = &ctx->forcing.UGRD_10maboveground;       return BMI_SUCCESS; }
    if (strcmp(name, "VGRD_10maboveground") == 0){*dest = &ctx->forcing.VGRD_10maboveground;       return BMI_SUCCESS; }
    if (strcmp(name, "param_catchment_vegetated_fraction") == 0) { *dest = &ctx->parameters.catchment_vegetated_fraction; return BMI_SUCCESS; }
    if (strcmp(name, "bare_soil_rsurf_exp") == 0)                { *dest = &ctx->parameters.bare_soil_rsurf_exp;          return BMI_SUCCESS; }
    if (strcmp(name, "verbosity") == 0)         { *dest = &ctx->options.verbosity;                  return BMI_SUCCESS; }
    /* --- calibration parameters --- */
    double *pp = param_field_ptr(ctx, name);
    if(pp != NULL) { 
        *dest = pp; 
        return BMI_SUCCESS; 
    }

    return BMI_FAILURE;
}

static int Get_value_at_indices(Bmi *self, const char *name, void *dest, int *inds, int count) {
    if (strcmp(name, "state_soil_moisture_theta") == 0 && CONTEXT(self)) {
        double *d = (double*)dest;
        for (int i = 0; i < count; i++) {
            if (inds[i] < 0 || inds[i] >= NDISC) return BMI_FAILURE;
            d[i] = CONTEXT(self)->state.soil_discrete_storage_theta[inds[i]];
        }
        return BMI_SUCCESS;
    }
    if (count == 1) return Get_value(self, name, dest);
    return BMI_FAILURE;
}

/* ================================================================== */
/*  Set value / Set value at indices                                   */
/* ================================================================== */

static int Set_value(Bmi *self, const char *name, void *src) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;

    /* Serialization protocol */
    if (strcmp(name, NGEN_SERIALIZATION_CREATE) == 0) {
        cfe_serialize_create(CONTEXT(self));
        return BMI_SUCCESS;
    }
    if (strcmp(name, NGEN_SERIALIZATION_FREE) == 0) {
        cfe_serialize_free(CONTEXT(self));
        return BMI_SUCCESS;
    }
    if (strcmp(name, NGEN_SERIALIZATION_STATE) == 0) {
        return cfe_serialize_deserialize(CONTEXT(self), (const char*)src);
    }

    /* Arrays need element-by-element copy into the target array */
    if (strcmp(name, "state_soil_moisture_theta") == 0) {
        double *s = (double*)src;
        for (int i = 0; i < NDISC; i++)
            CONTEXT(self)->state.soil_discrete_storage_theta[i] = s[i];
        return BMI_SUCCESS;
    }
    if (strcmp(name, "state_nash_subsurface_storage") == 0) {
        double *s = (double*)src;
        for (int i = 0; i < 2; i++)
            CONTEXT(self)->state.nash_subsurface_storage_m[i] = s[i];
        return BMI_SUCCESS;
    }
    if (strcmp(name, "state_giuh_queue") == 0) {
        double *s = (double*)src;
        for (int i = 0; i < CONTEXT(self)->parameters.giuh_num_ordinates; i++)
            CONTEXT(self)->state.giuh_queue_m[i] = s[i];
        return BMI_SUCCESS;
    }
    /* Params with BMI-boundary unit conversion (user-facing → internal SI) */
    if (strcmp(name, "soil_saturated_hydraulic_conductivity") == 0) {
        CONTEXT(self)->parameters.ksat_m_per_s = cm_per_h_to_m_per_s(*(double*)src);
        CONTEXT(self)->params_dirty = 1;
        return BMI_SUCCESS;
    }
    if (strcmp(name, "soil_saturated_capillary_head") == 0) {
        CONTEXT(self)->parameters.sat_capillary_head_m = cm_to_m(*(double*)src);
        CONTEXT(self)->params_dirty = 1;
        return BMI_SUCCESS;
    }
    /* Vegetated fraction: must re-derive bare_soil_fraction */
    if (strcmp(name, "param_catchment_vegetated_fraction") == 0) {
        double value = *(double*)src;
        if (value < 0.0 || value > 1.0) return BMI_FAILURE;
        CONTEXT(self)->parameters.catchment_vegetated_fraction = value;
        return (update_catchment_land_cover_fractions(&CONTEXT(self)->parameters) == 0)
               ? BMI_SUCCESS : BMI_FAILURE;
    }
    /* Bare-soil surface resistance exponent (Sakaguchi & Zeng) */
    if (strcmp(name, "bare_soil_rsurf_exp") == 0) {
        double value = *(double*)src;
        if (value <= 0.0) return BMI_FAILURE;
        CONTEXT(self)->parameters.bare_soil_rsurf_exp = value;
        return BMI_SUCCESS;
    }

    /* day_of_year from an external module suppresses internal calculation */
    if (strcmp(name, "day_of_year") == 0) {
        CONTEXT(self)->forcing.day_of_year = *(int*)src;
        CONTEXT(self)->forcing.day_of_year_set_externally = TRUE;
        return BMI_SUCCESS;
    }

    /* All other scalars: delegate through get_value_ptr */
    void *ptr = NULL;
    if (Get_value_ptr(self, name, &ptr) != BMI_SUCCESS || ptr == NULL)
        return BMI_FAILURE;

    int nbytes = 0;
    if (Get_var_nbytes(self, name, &nbytes) != BMI_SUCCESS)
        return BMI_FAILURE;

    memcpy(ptr, src, nbytes);

    /* Mark dirty if this was a calibration parameter */
    if (param_field_ptr(CONTEXT(self), name) != NULL)
        CONTEXT(self)->params_dirty = 1;

    return BMI_SUCCESS;
}

static int Set_value_at_indices(Bmi *self, const char *name, int *inds, int count, void *src) {
    if (strcmp(name, "state_soil_moisture_theta") == 0 && CONTEXT(self)) {
        double *s = (double*)src;
        for (int i = 0; i < count; i++) {
            if (inds[i] < 0 || inds[i] >= NDISC) return BMI_FAILURE;
            CONTEXT(self)->state.soil_discrete_storage_theta[inds[i]] = s[i];
        }
        return BMI_SUCCESS;
    }
    if (count == 1) return Set_value(self, name, src);
    return BMI_FAILURE;
}

/* ================================================================== */
/*  Grid information                                                   */
/* ================================================================== */

static int Get_grid_rank(Bmi *self, int grid, int *rank) {
    *rank = (grid == 0) ? 0 : 1;
    return BMI_SUCCESS;
}

static int Get_grid_size(Bmi *self, int grid, int *size) {
    CFE_Model_Context *ctx = CONTEXT(self);
    switch (grid) {
        case 0: *size = 1; break;
        case 1: *size = NDISC; break;
        case 2: *size = 2; break;  /* subsurface Nash (always 2) */
        case 3:
            if (ctx)
                *size = ctx->parameters.giuh_num_ordinates;
            else
                *size = 0;
            break;
        default: return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

static int Get_grid_type(Bmi *self, int grid, char *type) {
    if (grid < 0 || grid > 3) return BMI_FAILURE;
    strcpy(type, (grid == 0) ? "scalar" : "vector");
    return BMI_SUCCESS;
}

static int Get_grid_shape(Bmi *self, int grid, int *shape) {
    if (grid > 0) return Get_grid_size(self, grid, shape);
    return BMI_FAILURE;
}

static int Get_grid_spacing  (Bmi *self, int g, double *s)       { return BMI_FAILURE; }
static int Get_grid_origin   (Bmi *self, int g, double *o)       { return BMI_FAILURE; }
static int Get_grid_x        (Bmi *self, int g, double *x)       { return BMI_FAILURE; }
static int Get_grid_y        (Bmi *self, int g, double *y)       { return BMI_FAILURE; }
static int Get_grid_z        (Bmi *self, int g, double *z)       { return BMI_FAILURE; }
static int Get_grid_node_count(Bmi *self, int g, int *c)         { return Get_grid_size(self, g, c); }
static int Get_grid_edge_count(Bmi *self, int g, int *c)         { return BMI_FAILURE; }
static int Get_grid_face_count(Bmi *self, int g, int *c)         { return BMI_FAILURE; }
static int Get_grid_edge_nodes(Bmi *self, int g, int *e)         { return BMI_FAILURE; }
static int Get_grid_face_edges(Bmi *self, int g, int *e)         { return BMI_FAILURE; }
static int Get_grid_face_nodes(Bmi *self, int g, int *n)         { return BMI_FAILURE; }
static int Get_grid_nodes_per_face(Bmi *self, int g, int *n)     { return BMI_FAILURE; }

/* ================================================================== */
/*  Registration                                                       */
/* ================================================================== */

Bmi* register_bmi_cfe(Bmi *model) {
    if (model == NULL) return NULL;

    model->data = NULL;

    model->initialize          = Initialize;
    model->update              = Update;
    model->update_until        = Update_until;
    model->finalize            = Finalize;

    model->get_component_name  = Get_component_name;
    model->get_input_item_count  = Get_input_item_count;
    model->get_output_item_count = Get_output_item_count;
    model->get_input_var_names   = Get_input_var_names;
    model->get_output_var_names  = Get_output_var_names;

    model->get_var_grid        = Get_var_grid;
    model->get_var_type        = Get_var_type;
    model->get_var_units       = Get_var_units;
    model->get_var_itemsize    = Get_var_itemsize;
    model->get_var_nbytes      = Get_var_nbytes;
    model->get_var_location    = Get_var_location;

    model->get_current_time    = Get_current_time;
    model->get_start_time      = Get_start_time;
    model->get_end_time        = Get_end_time;
    model->get_time_units      = Get_time_units;
    model->get_time_step       = Get_time_step;

    model->get_value           = Get_value;
    model->get_value_ptr       = Get_value_ptr;
    model->get_value_at_indices = Get_value_at_indices;

    model->set_value           = Set_value;
    model->set_value_at_indices = Set_value_at_indices;

    model->get_grid_rank       = Get_grid_rank;
    model->get_grid_size       = Get_grid_size;
    model->get_grid_type       = Get_grid_type;
    model->get_grid_shape      = Get_grid_shape;
    model->get_grid_spacing    = Get_grid_spacing;
    model->get_grid_origin     = Get_grid_origin;
    model->get_grid_x          = Get_grid_x;
    model->get_grid_y          = Get_grid_y;
    model->get_grid_z          = Get_grid_z;
    model->get_grid_node_count = Get_grid_node_count;
    model->get_grid_edge_count = Get_grid_edge_count;
    model->get_grid_face_count = Get_grid_face_count;
    model->get_grid_edge_nodes = Get_grid_edge_nodes;
    model->get_grid_face_edges = Get_grid_face_edges;
    model->get_grid_face_nodes = Get_grid_face_nodes;
    model->get_grid_nodes_per_face = Get_grid_nodes_per_face;

    return model;
}

CFE_Model_Context* new_bmi_cfe(void) {
    /* Allocate a zeroed context — Initialize will populate it */
    return (CFE_Model_Context*)calloc(1, sizeof(CFE_Model_Context));
}

void delete_bmi_cfe(Bmi *model) {
    if (model != NULL) {
        if (model->data != NULL) {
            Finalize(model);
        }
        free(model);
    }
}
