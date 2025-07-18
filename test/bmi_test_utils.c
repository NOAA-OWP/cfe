#include "bmi_test_utils.h"
#include "general_test_utils.h"

TestFixture* setup(const char* cfg_file)
{
    TestFixture* fixture = malloc(sizeof(TestFixture));
    fixture->bmi_model = (Bmi *) malloc(sizeof(Bmi));
    fixture->cfg_file = cfg_file;

    char* var_names[EXPECTED_TOTAL_VAR_COUNT] = {
        "RAIN_RATE",
        "GIUH_RUNOFF",
        "INFILTRATION_EXCESS",
        "DIRECT_RUNOFF",
        "NASH_LATERAL_RUNOFF",
        "DEEP_GW_TO_CHANNEL_FLUX",
        "SOIL_TO_GW_FLUX",
        "Q_OUT",
        "POTENTIAL_ET",
        "ACTUAL_ET",
        "GW_STORAGE",
        "SOIL_STORAGE",
        "SOIL_STORAGE_CHANGE",
        "SURF_RUNOFF_SCHEME",
        "NWM_PONDED_DEPTH",
        "atmosphere_water__liquid_equivalent_precipitation_rate",
        "water_potential_evaporation_flux",
        "ice_fraction_schaake",
        "ice_fraction_xinanjiang",
        "soil_moisture_profile"
    };

    fixture->expected_output_and_input_var_names = allocate_array_of_strings(EXPECTED_TOTAL_VAR_COUNT, BMI_MAX_VAR_NAME);
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        strcpy(fixture->expected_output_and_input_var_names[i], var_names[i]);
    }

    fixture->expected_output_var_names = fixture->expected_output_and_input_var_names;
    fixture->expected_input_var_names = fixture->expected_output_and_input_var_names + EXPECTED_OUTPUT_VAR_COUNT;

    // Note that for now, grid id for all variables is 0
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        fixture->expected_grid_ids[i] = 0;
    }

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
void get_arbitrary_input_var_values(int example_case, double current_model_time, double* value_array) {
    // For now, use the same simple group of values for everything
    // TODO: might need to confirm the validity (or the ideal-ness) of these values further
    double arbitrary_input_var_values[EXPECTED_INPUT_VAR_COUNT] = {0.55, 0.27, 0.1, 0.1, 0.2};
    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++)
        value_array[i] = arbitrary_input_var_values[i];
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
 * @param example_case The specific example test case, which could affect which values are used.
 * @param current_model_time The current model time, which could affect which values are used.
 * @return Whether the set operation was successful.
 */
bool set_module_input_variables_before_update(const TestFixture* fixture, const int example_case, const double current_model_time) {
    int bmi_status;
    char var_type[BMI_MAX_TYPE_NAME];

    double arbitrary_input_var_values[EXPECTED_INPUT_VAR_COUNT];
    get_arbitrary_input_var_values(example_case, current_model_time, arbitrary_input_var_values);

    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++) {
        // Sanity check
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, fixture->expected_input_var_names[i], var_type);
        if (bmi_status == BMI_FAILURE) {
            printf("\nCan't set module inputs to advance; test helper function encountered BMI_FAILURE getting type of variable '%s' for sanity check", fixture->expected_input_var_names[i]);
            return false;
        }
        if (!confirm_matches_expected_strs("double", var_type)) {
            printf("\nCan't set module inputs to advance; test helper function has type mismatch for variable %s", fixture->expected_input_var_names[i]);
            return false;
        }
        // Assuming the sanity check of the type is good, set the above-prepared arbitrary value for this
        bmi_status = fixture->bmi_model->set_value(fixture->bmi_model, fixture->expected_input_var_names[i], arbitrary_input_var_values + i);
        if (bmi_status == BMI_FAILURE) {
            printf("\nCan't set module inputs to advance; test helper function encountered BMI_FAILURE attempting to set variable '%s'", fixture->expected_input_var_names[i]);
            return false;
        }
    }
    return true;
}