#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmi_test_utils.h"
#include "general_test_utils.h"

TestFixture* setup(const unsigned int example_case, const char* cfg_file)
{
    TestFixture* fixture = malloc(sizeof(TestFixture));
    fixture->bmi_model = (Bmi *) malloc(sizeof(Bmi));
    fixture->current_test_example = example_case;
    fixture->cfg_file = cfg_file;

    char* var_names[EXPECTED_TOTAL_VAR_COUNT] = {
        /* 32 outputs */
        "discharge_m",
        "surface_runoff_m",
        "lateral_flow_m",
        "baseflow_m",
        "actual_et_m",
        "vol_balance_residual_m",
        "state_soil_storage_m",
        "state_gw_storage_m",
        "state_current_timestep",
        "state_soil_moisture_theta",
        "state_nash_subsurface_storage",
        "state_giuh_queue",
        "config_simulate_discrete_soil_moisture",
        "param_catchment_area_km2",
        "param_soil_depth_m",
        "param_soil_porosity",
        "timestep_storage_start_m",
        "timestep_input_m",
        "timestep_output_m",
        "timestep_storage_end_m",
        "potential_et_m",
        "giuh_outflow_m",
        "soil_to_gw_percolation_flux_m",
        "soil_moisture_theta_1",
        "soil_moisture_theta_2",
        "soil_moisture_theta_3",
        "soil_moisture_theta_4",
        "bare_soil_evaporation_m",
        "impervious_runoff_m",
        "pervious_runoff_m",
        "surface_routed_to_outlet_m",
        "lateral_flow_generated_m",
        /* 5 inputs (PT and soil evap disabled in test config; AORC fields omitted) */
        "rainfall_depth_m",
        "et_potential_m",
        "ice_fraction",
        "param_catchment_vegetated_fraction",
        "bare_soil_rsurf_exp"
    };

    fixture->expected_output_and_input_var_names = allocate_array_of_strings(EXPECTED_TOTAL_VAR_COUNT, BMI_MAX_VAR_NAME);
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        strcpy(fixture->expected_output_and_input_var_names[i], var_names[i]);
    }

    fixture->expected_output_var_names = fixture->expected_output_and_input_var_names;
    fixture->expected_input_var_names = fixture->expected_output_and_input_var_names + EXPECTED_OUTPUT_VAR_COUNT;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++)
        fixture->expected_grid_ids[i] = 0;
    /* array outputs have non-zero grid ids */
    fixture->expected_grid_ids[9]  = 1;  /* state_soil_moisture_theta */
    fixture->expected_grid_ids[10] = 2;  /* state_nash_subsurface_storage */
    fixture->expected_grid_ids[11] = 3;  /* state_giuh_queue */

    register_bmi_cfe(fixture->bmi_model);

    return fixture;
}

void teardown(TestFixture* fixture)
{
    free_array_of_strings(fixture->expected_output_and_input_var_names, EXPECTED_TOTAL_VAR_COUNT);
    free(fixture->bmi_model);
}

/**
 * Helper to dynamically allocate memory and get all BMI variable names, plus the variable counts.
 *
 * Regarding the ordering of variable names in the returned array, it should be the continuous collection of all output
 * variables followed by all input variables.  Within each of those continuous collections, the variables should be
 * ordered in the same ways as they are when returned by ``get_output_var_names`` and ``get_input_var_names``.
 *
 * Note that if the function succeeds, it will allocate memory in ``all_var_names`` that must later be freed.  However,
 * if it does not succeed (e.g., the attempt to get the output variable names from the model fails), even if it
 * allocated memory momentarily, it will free that memory before returning.
 *
 * @param bmi_model Pointer to the BMI model itself
 * @param input_var_count Allocated memory for the number of output variables
 * @param output_var_count Allocated memory for the number input variables
 * @return Allocate pointer(s) containing an array of strings with variable names (output variables first), or NULL on failure
 */
char** get_all_bmi_variable_names(Bmi* bmi_model, int* output_var_count, int* input_var_count)
{
    int bmi_status = bmi_model->get_output_item_count(bmi_model, output_var_count);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status attempting to get output variable count (in order to get all names)");
        return NULL;
    }
    bmi_status = bmi_model->get_input_item_count(bmi_model, input_var_count);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status attempting to get input variable count (in order to get all names)");
        return NULL;
    }
    int total_var_count = *input_var_count + *output_var_count;

    // Sanity check that there is something
    if (total_var_count == 0) {
        printf("\nWhile calls to get output/input variable counts succeeded, total variable count was 0");
        return NULL;
    }

    // *** IMPORTANT *** - Now that this is done, any failure return must be proceeded by freeing this memory
    char** names = allocate_array_of_strings(total_var_count, BMI_MAX_VAR_NAME);

    bmi_status = bmi_model->get_output_var_names(bmi_model, names);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code getting output variable names (while getting all names");
        free_array_of_strings(names, total_var_count);
        return NULL;
    }

    // Do some pointer arithmatic here to start after the output names in all_var_names
    bmi_status = bmi_model->get_input_var_names(bmi_model, names + *output_var_count);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code getting input variable names (while getting all names)");
        free_array_of_strings(names, total_var_count);
        return NULL;
    }

    return names;
}

/**
 * Get an array of arbitrary but valid values to use to set for module inputs, and save them to a provided array.
 *
 * @param example_case The specific example test case, which could affect which values are used.
 * @param current_model_time The current model time, which could affect which values are used.
 * @param value_array Pointer to the array in which to save the values (which must be of size EXPECTED_INPUT_VAR_COUNT).
 */
void get_arbitrary_input_var_values(const unsigned int example_case, double current_model_time, double* value_array) {
    // For now, use the same simple group of values for everything
    // TODO: might need to confirm the validity (or the ideal-ness) of these values further
    /* v3 with AORC disabled: rainfall_depth_m, et_potential_m, ice_fraction,
     *     param_catchment_vegetated_fraction, bare_soil_rsurf_exp */
    double arbitrary_input_var_values[EXPECTED_INPUT_VAR_COUNT] = {
        0.001, 0.0001, 0.0,
        0.7, 5.0
    };
    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++)
        value_array[i] = arbitrary_input_var_values[i];
}

/**
 * Get values of output variables, capturing them in a provided double array (casting when needed).
 *
 * @param fixture The test fixture, with the module and array of output variable names.
 * @param value_array The array of doubles in which to place values.
 * @return Whether the operation executed successfully.
 */
bool get_output_var_values(TestFixture* fixture, double* value_array)
{
    char var_type[BMI_MAX_VAR_NAME];
    int bmi_status, int_var_val;
    double double_var_val;

    for (int i = 0; i < EXPECTED_OUTPUT_VAR_COUNT; i++) {
        // Have local var for these just for readability
        const char* var_name = fixture->expected_output_var_names[i];
        /* Determine type dynamically from the model */
        char check_type[BMI_MAX_TYPE_NAME];
        fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, check_type);
        bool is_int = (strcmp(check_type, "int") == 0);
        bool is_double = (strcmp(check_type, "double") == 0);

        /* skip array variables (grid id > 0) and non-numeric types */
        int grid_id = 0;
        fixture->bmi_model->get_var_grid(fixture->bmi_model, var_name, &grid_id);
        if (grid_id > 0) {
            printf("\n  [get_output_var_values] skipping array variable '%s' (grid %d)", var_name, grid_id);
            value_array[i] = 0.0;
            continue;
        }
        if (!is_int && !is_double) {
            printf("\n  [get_output_var_values] skipping non-numeric variable '%s' (type '%s')", var_name, check_type);
            value_array[i] = 0.0;
            continue;
        }

        void* val_ptr = is_double ? (void*) &double_var_val : (void*) &int_var_val;
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, val_ptr);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting uninitialized value for output '%s'", var_name);
            return false;
        }
        value_array[i] = is_int ? (double)int_var_val : double_var_val;
    }
    return true;
}

/*
 * Setup expected values for grid ids for all the output and input (in that order) BMI variables.
 *
 * The array is assumed to be of a hardcoded size (``EXPECTED_TOTAL_VAR_COUNT``).
 *
 * Regarding ordering or the setup array, the value set at any index `n` should be the grid id for the variable at index
 * `n` in the ``get_all_bmi_variable_names`` function.
 *
 * @param grid_id_array Address to start of array in which to save grid ids, going in the same order as variable names
 *                      returned by BMI getters, with all output variable names first followed by input names.
 */
//void setup_expected_grid_ids(int* grid_id_array)
//{
//    // Note that for now, grid id for all variables is 0
//    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
//        grid_id_array[i] = 0;
//    }
//}

/**
 * Set all necessary module BMI input variables to reasonable values, as needed prior to advancing the model.
 *
 * @param fixture The test fixture, which contains the module.
 * @param current_model_time The current model time, which could affect which values are used.
 * @return Whether the set operation was successful.
 */
bool set_arbitrary_input_variables_before_update(const TestFixture* fixture, const double current_model_time)
{
    double arbitrary_input_var_values[EXPECTED_INPUT_VAR_COUNT];
    get_arbitrary_input_var_values(fixture->current_test_example, current_model_time, arbitrary_input_var_values);
    return set_specified_input_variables_before_update(fixture, current_model_time, arbitrary_input_var_values);
}

/**
 * Set all necessary module BMI input variables to specified values, as needed prior to advancing the model.
 *
 * @param fixture The test fixture, which contains the module.
 * @param current_model_time The current model time, which could affect which values are used
 * @param input_var_values The values to use to set BMI input variables, ordered in the same way as the variable names
 * when retrieved.
 * @return Whether the set operation was successful.
 */
bool set_specified_input_variables_before_update(const TestFixture* fixture, double current_model_time, double* input_var_values) {
    int bmi_status;
    char var_type[BMI_MAX_TYPE_NAME];

    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++) {
        // Sanity check
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, fixture->expected_input_var_names[i], var_type);
        if (bmi_status == BMI_FAILURE) {
            printf("\nCan't set module inputs to advance; test helper function encountered BMI_FAILURE getting type of variable '%s' for sanity check", fixture->expected_input_var_names[i]);
            return false;
        }
        /* v3: only set double inputs; log skips for non-double types */
        if (strcmp(var_type, "double") != 0) {
            printf("\n  [set_inputs] skipping non-double input '%s' (type '%s')",
                   fixture->expected_input_var_names[i], var_type);
            continue;
        }
        bmi_status = fixture->bmi_model->set_value(fixture->bmi_model, fixture->expected_input_var_names[i], input_var_values + i);
        if (bmi_status == BMI_FAILURE) {
            printf("\nCan't set module inputs to advance; test helper function encountered BMI_FAILURE attempting to set variable '%s'", fixture->expected_input_var_names[i]);
            return false;
        }
    }
    return true;
}