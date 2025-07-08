#include "general_test_utils.h"
#include "bmi_test_utils.h"
#include "bmi_cfe.h"


int test_get_component_name(TestFixture* fixture)
{
    char name[BMI_MAX_COMPONENT_NAME];
    int bmi_status = fixture->bmi_model->get_component_name(fixture->bmi_model, name);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code");
        return TEST_RETURN_CODE_FAIL;
    }
    if (confirm_matches_expected_strs(EXPECTED_COMPONENT_NAME, name))
        return TEST_RETURN_CODE_PASS;
    return TEST_RETURN_CODE_FAIL;
}

int test_get_grid_rank(TestFixture* fixture)
{
    // Can rely on expected id values and assume if they are off it's caught by test_get_var_grid
    int ordered_grid_ids[EXPECTED_TOTAL_VAR_COUNT];
    // This assumes the passed array will be of size EXPECTED_TOTAL_VAR_COUNT
    setup_expected_grid_ids(&ordered_grid_ids);

    // However, still set this up for more clear messaging when there are errors
    int output_count, input_count, actual_rank, status = 1;

    char** all_var_names = get_all_bmi_variable_names(fixture->bmi_model, &output_count, &input_count);
    if (all_var_names == NULL) {
        printf("\nFailed getting variable names to test get_grid_rank");
        return TEST_RETURN_CODE_FAIL;
    }
    // *** IMPORTANT *** - IFF get_all_bmi_variable_names was successful (i.e., we get this far), memory will be
    // allocated in all_var_names that later needs to be freed (via free_array_of_strings)
    int total_var_count = output_count + input_count;

    // Sanity check this here
    if (total_var_count != EXPECTED_TOTAL_VAR_COUNT) {
        printf("\nUnexpected total variable count (%i, but expected %i) getting variable names to test get_grid_rank",
               total_var_count, EXPECTED_TOTAL_VAR_COUNT);
        free_array_of_strings(all_var_names, total_var_count);
        return TEST_RETURN_CODE_FAIL;
    }

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_grid_rank(fixture->bmi_model, ordered_grid_ids[i], &actual_rank);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting grid rank for '%s'", all_var_names[i]);
            free_array_of_strings(all_var_names, total_var_count);
            return TEST_RETURN_CODE_FAIL;
        }

        // For now, all grids are of rank 1
        if (!confirm_matches_expected_ints(1, actual_rank)) {
            printf("\nGrid rank for '%s' did not match expected", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }
        status = TEST_RETURN_CODE_PASS;
    }
    free_array_of_strings(all_var_names, total_var_count);
    return status;



}

int test_get_input_item_count(TestFixture* fixture)
{
    int input_count, bmi_status;
    bmi_status = fixture->bmi_model->get_input_item_count(fixture->bmi_model, &input_count);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code");
        return TEST_RETURN_CODE_FAIL;
    }
    if (confirm_matches_expected_ints(EXPECTED_INPUT_VAR_COUNT, input_count))
        return TEST_RETURN_CODE_PASS;
    return TEST_RETURN_CODE_FAIL;
}

int test_get_input_var_names(TestFixture* fixture)
{
    int var_count, bmi_status, result;
    bmi_status = fixture->bmi_model->get_input_item_count(fixture->bmi_model, &var_count);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get input variable count (in order to get names)");
        return TEST_RETURN_CODE_FAIL;
    }

    char** actual_var_names = allocate_array_of_strings(var_count, BMI_MAX_VAR_NAME);

    const char* expected_var_names[EXPECTED_INPUT_VAR_COUNT] = {
        "atmosphere_water__liquid_equivalent_precipitation_rate",
        "water_potential_evaporation_flux",
        "ice_fraction_schaake",
        "ice_fraction_xinanjiang",
        "soil_moisture_profile"
    };

    bmi_status = fixture->bmi_model->get_input_var_names(fixture->bmi_model, actual_var_names);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code getting input variable names");
        result = TEST_RETURN_CODE_FAIL;
    }
    else if (confirm_matches_expected_str_arrays(expected_var_names, actual_var_names, var_count)) {
        result = TEST_RETURN_CODE_PASS;
    }
    else {
        result = TEST_RETURN_CODE_FAIL;
    }
    free_array_of_strings(actual_var_names, var_count);
    return result;
}

int test_get_output_item_count(TestFixture* fixture)
{
    int output_count, bmi_status;
    bmi_status = fixture->bmi_model->get_output_item_count(fixture->bmi_model, &output_count);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code");
        return TEST_RETURN_CODE_FAIL;
    }
    if (confirm_matches_expected_ints(EXPECTED_OUTPUT_VAR_COUNT, output_count))
        return TEST_RETURN_CODE_PASS;
    return TEST_RETURN_CODE_FAIL;
}

int test_get_output_var_names(TestFixture* fixture)
{
    int var_count, bmi_status, result;
    bmi_status = fixture->bmi_model->get_output_item_count(fixture->bmi_model, &var_count);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get output variable count (in order to get names)");
        return TEST_RETURN_CODE_FAIL;
    }

    char** actual_var_names = allocate_array_of_strings(var_count, BMI_MAX_VAR_NAME);

    const char* expected_var_names[EXPECTED_OUTPUT_VAR_COUNT] = {
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
        "NWM_PONDED_DEPTH"
    };

    bmi_status = fixture->bmi_model->get_output_var_names(fixture->bmi_model, actual_var_names);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code getting output variable names");
        result = TEST_RETURN_CODE_FAIL;
    }
    else if (confirm_matches_expected_str_arrays(expected_var_names, actual_var_names, var_count)) {
        result = TEST_RETURN_CODE_PASS;
    }
    else {
        result = TEST_RETURN_CODE_FAIL;
    }
    free_array_of_strings(actual_var_names, var_count);
    return result;
}

int test_get_var_grid(TestFixture* fixture)
{
    int output_count, input_count, grid_value, status;
    char** all_var_names = get_all_bmi_variable_names(fixture->bmi_model, &output_count, &input_count);
    if (all_var_names == NULL) {
        printf("\nFailed getting variable names to test get_var_grid");
        return TEST_RETURN_CODE_FAIL;
    }
    // *** IMPORTANT *** - IFF get_all_bmi_variable_names was successful (i.e., we get this far), memory will be
    // allocated in all_var_names that later needs to be freed
    int total_var_count = output_count + input_count;

    // Sanity check this here
    if (total_var_count != EXPECTED_TOTAL_VAR_COUNT) {
        printf("\nUnexpected total variable count (%i, but expected %i) getting variable names to test get_var_grid",
               total_var_count, EXPECTED_TOTAL_VAR_COUNT);
        return TEST_RETURN_CODE_FAIL;
    }

    // Note that all grid values are 0, but we'll do this to stick to the pattern
    int ordered_grid_expected_values[EXPECTED_TOTAL_VAR_COUNT];
    // This assumes the passed array will be of size EXPECTED_TOTAL_VAR_COUNT
    setup_expected_grid_ids(&ordered_grid_expected_values);

    for (int i = 0; i < total_var_count; i++) {
        status = fixture->bmi_model->get_var_grid(fixture->bmi_model, all_var_names[i], &grid_value);
        if (status == BMI_FAILURE) {
            printf("\nReturned BMI_FAILURE status code getting grid for '%s'", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }
        if (!confirm_matches_expected_ints(ordered_grid_expected_values[i], grid_value)) {
            status = TEST_RETURN_CODE_FAIL;
            printf("\nGrid value for are different for %s", all_var_names[i]);
            break;
        }
        status = TEST_RETURN_CODE_PASS;
    }
    free_array_of_strings(all_var_names, total_var_count);
    return status;
}

int test_get_var_itemsize(TestFixture* fixture)
{
    int output_count, input_count, item_size, expected_size, status = 1;
    char** all_var_names = get_all_bmi_variable_names(fixture->bmi_model, &output_count, &input_count);
    if (all_var_names == NULL) {
        printf("\nFailed getting variable names to test get_var_itemsize");
        return TEST_RETURN_CODE_FAIL;
    }
    // *** IMPORTANT *** - IFF get_all_bmi_variable_names was successful (i.e., we get this far), memory will be
    // allocated in all_var_names that later needs to be freed
    int total_var_count = output_count + input_count;

    for (int i = 0; i < total_var_count; i++) {
        status = fixture->bmi_model->get_var_itemsize(fixture->bmi_model, all_var_names[i], &item_size);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting item size for '%s'", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }

        // Only SURF_RUNOFF_SCHEME is int, and the rest (output and input) are doubles
        expected_size = strcmp(all_var_names[i], "SURF_RUNOFF_SCHEME") == 0 ? sizeof(int) : sizeof(double);

        if (!confirm_matches_expected_ints(expected_size, item_size)) {
            printf("\nSize for '%s' did not match expected", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }
        // Otherwise ...
        status = TEST_RETURN_CODE_PASS;
    }
    free_array_of_strings(all_var_names, total_var_count);
    return status;
}

int test_get_var_location(TestFixture* fixture)
{
    int output_count, input_count, status = 1;
    char actual_value[BMI_MAX_VAR_NAME];
    char** all_var_names = get_all_bmi_variable_names(fixture->bmi_model, &output_count, &input_count);
    if (all_var_names == NULL) {
        printf("\nFailed getting variable names to test get_var_location");
        return TEST_RETURN_CODE_FAIL;
    }
    // *** IMPORTANT *** - IFF get_all_bmi_variable_names was successful (i.e., we get this far), memory will be
    // allocated in all_var_names that later needs to be freed
    int total_var_count = output_count + input_count;

    for (int i = 0; i < total_var_count; i++) {
        status = fixture->bmi_model->get_var_location(fixture->bmi_model, all_var_names[i], &actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting location for '%s'", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }

        // Only SURF_RUNOFF_SCHEME is "none", and the rest (output and input) are "node"
        bool is_runoff_scheme_var = strcmp(all_var_names[i], "SURF_RUNOFF_SCHEME") == 0;

        if (is_runoff_scheme_var && confirm_matches_expected_strs("none", actual_value))
            status = TEST_RETURN_CODE_PASS;
        else if (!is_runoff_scheme_var && confirm_matches_expected_strs("node", actual_value))
            status = TEST_RETURN_CODE_PASS;
        else
            status = TEST_RETURN_CODE_FAIL;

        if (status != TEST_RETURN_CODE_PASS) {
            printf("\nLocation for '%s' did not match expected", all_var_names[i]);
            break;
        }
    }

    free_array_of_strings(all_var_names, total_var_count);
    return status;
}


int test_get_var_units(TestFixture* fixture)
{
    int output_count, input_count, status = 1;
    char actual_value[BMI_MAX_VAR_NAME];
    char** all_var_names = get_all_bmi_variable_names(fixture->bmi_model, &output_count, &input_count);
    if (all_var_names == NULL) {
        printf("\nFailed getting variable names to test get_var_units");
        return TEST_RETURN_CODE_FAIL;
    }
    // *** IMPORTANT *** - IFF get_all_bmi_variable_names was successful (i.e., we get this far), memory will be
    // allocated in all_var_names that later needs to be freed
    int total_var_count = output_count + input_count;

    char* expected_unit_values[EXPECTED_TOTAL_VAR_COUNT] = {
        "m",        // RAIN_RATE
        "m",        // GIUH_RUNOFF
        "m",        // INFILTRATION_EXCESS
        "m",        // DIRECT_RUNOFF
        "m",        // NASH_LATERAL_RUNOFF
        "m",        // DEEP_GW_TO_CHANNEL_FLUX
        "m",        // SOIL_TO_GW_FLUX
        "m",        // Q_OUT
        "m",        // POTENTIAL_ET
        "m",        // ACTUAL_ET
        "m",        // GW_STORAGE
        "m",        // SOIL_STORAGE
        "m",        // SOIL_STORAGE_CHANGE
        "none",     // SURF_RUNOFF_SCHEME
        "m",        // NWM_PONDED_DEPTH
        "mm h-1",   // atmosphere_water__liquid_equivalent_precipitation_rate
        "m s-1",    // water_potential_evaporation_flux
        "m",        // ice fraction in meters
        "none",     // ice fraction [-]
        "none"      // soil moisture profile is in decimal fraction -rlm
    };

    for (int i = 0; i < total_var_count; i++) {
        status = fixture->bmi_model->get_var_units(fixture->bmi_model, all_var_names[i], &actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting units for '%s'", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }
        if (!confirm_matches_expected_strs(expected_unit_values[i], actual_value)) {
            printf("\nUnits for '%s' did not match expected", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }
        status = TEST_RETURN_CODE_PASS;
    }

    free_array_of_strings(all_var_names, total_var_count);
    return status;
}

int test_get_var_type(TestFixture* fixture)
{
    int output_count, input_count, status = 1;
    char actual_value[BMI_MAX_TYPE_NAME];
    char** all_var_names = get_all_bmi_variable_names(fixture->bmi_model, &output_count, &input_count);
    if (all_var_names == NULL) {
        printf("\nFailed getting variable names to test get_var_type");
        return TEST_RETURN_CODE_FAIL;
    }
    // *** IMPORTANT *** - IFF get_all_bmi_variable_names was successful (i.e., we get this far), memory will be
    // allocated in all_var_names that later needs to be freed
    int total_var_count = output_count + input_count;

    for (int i = 0; i < total_var_count; i++) {
        status = fixture->bmi_model->get_var_type(fixture->bmi_model, all_var_names[i], &actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting type for '%s'", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }

        // Only SURF_RUNOFF_SCHEME is int, and the rest (output and input) are doubles
        bool is_runoff_scheme_var = strcmp(all_var_names[i], "SURF_RUNOFF_SCHEME") == 0;

        if (is_runoff_scheme_var && confirm_matches_expected_strs("int", actual_value))
            status = TEST_RETURN_CODE_PASS;
        else if (!is_runoff_scheme_var && confirm_matches_expected_strs("double", actual_value))
            status = TEST_RETURN_CODE_PASS;
        else {
            printf("\nType for '%s' did not match expected", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }
    }

    free_array_of_strings(all_var_names, total_var_count);
    return status;
}

int test_get_var_nbytes(TestFixture* fixture)
{
    int output_count, input_count, item_nbytes, expected_nbytes, status = 1;
    char** all_var_names = get_all_bmi_variable_names(fixture->bmi_model, &output_count, &input_count);
    if (all_var_names == NULL) {
        printf("\nFailed getting variable names to test get_var_nbytes");
        return TEST_RETURN_CODE_FAIL;
    }
    // *** IMPORTANT *** - IFF get_all_bmi_variable_names was successful (i.e., we get this far), memory will be
    // allocated in all_var_names that later needs to be freed
    int total_var_count = output_count + input_count;

    for (int i = 0; i < total_var_count; i++) {
        status = fixture->bmi_model->get_var_nbytes(fixture->bmi_model, all_var_names[i], &item_nbytes);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting nbytes for '%s'", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }

        // Only SURF_RUNOFF_SCHEME is int, and the rest (output and input) are doubles
        // And everything should be of single element, so just the size
        expected_nbytes = strcmp(all_var_names[i], "SURF_RUNOFF_SCHEME") == 0 ? sizeof(int) : sizeof(double);

        if (!confirm_matches_expected_ints(expected_nbytes, item_nbytes)){
            printf("\nnbytes for '%s' did not match expected", all_var_names[i]);
            status = TEST_RETURN_CODE_FAIL;
            break;
        }
        status = TEST_RETURN_CODE_PASS;
    }
    free_array_of_strings(all_var_names, total_var_count);
    return status;
}

int test_initialize(TestFixture* fixture)
{
    if (fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file) == BMI_SUCCESS)
        return TEST_RETURN_CODE_PASS;
    return TEST_RETURN_CODE_FAIL;
}

int main(int argc, const char* argv[])
{
    char* config_file;
    int result = -1;

    // Test function is always argv[1]
    // Example test instance is argv[2] (or assumed to be "1" if not given)

    if (argc < 3 || strcmp(argv[2], "1") == 0) {
        config_file = BMI_INIT_CONFIG_EX_1;
    }
    else {
        printf("\nUnexpected test case %s\n", argv[2]);
        return 1;
    }

    TestFixture* fixture = setup(config_file);

    if (strcmp(argv[1], "test_initialize") == 0)
        result = test_initialize(fixture);
    // TODO: update test
    // TODO: update_until test
    // TODO: finalize test
    else if (strcmp(argv[1], "test_get_component_name") == 0)
        result = test_get_component_name(fixture);
    else if (strcmp(argv[1], "test_get_grid_rank") == 0)
        result = test_get_grid_rank(fixture);
    else if (strcmp(argv[1], "test_get_input_item_count") == 0)
        result = test_get_input_item_count(fixture);
    else if (strcmp(argv[1], "test_get_input_var_names") == 0)
        result = test_get_input_var_names(fixture);
    else if (strcmp(argv[1], "test_get_output_item_count") == 0)
        result = test_get_output_item_count(fixture);
    else if (strcmp(argv[1], "test_get_output_var_names") == 0)
        result = test_get_output_var_names(fixture);
    else if (strcmp(argv[1], "test_get_var_grid") == 0)
        result = test_get_var_grid(fixture);
    else if (strcmp(argv[1], "test_get_var_itemsize") == 0)
        result = test_get_var_itemsize(fixture);
    else if (strcmp(argv[1], "test_get_var_location") == 0)
        result = test_get_var_location(fixture);
    else if (strcmp(argv[1], "test_get_var_units") == 0)
        result = test_get_var_units(fixture);
    else if (strcmp(argv[1], "test_get_var_type") == 0)
        result = test_get_var_type(fixture);
    else if (strcmp(argv[1], "test_get_var_nbytes") == 0)
        result = test_get_var_nbytes(fixture);
    else
        printf("\nUnexpected test function %s\n", argv[1]);

    teardown(fixture);
    return result;
}