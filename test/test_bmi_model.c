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
    int status, actual_rank;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_grid_rank(fixture->bmi_model, fixture->expected_grid_ids[i], &actual_rank);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting grid rank for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // For now, all grids are of rank 1
        if (!confirm_matches_expected_ints(1, actual_rank)) {
            printf("\nGrid rank for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
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

    bmi_status = fixture->bmi_model->get_input_var_names(fixture->bmi_model, actual_var_names);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code getting input variable names");
        result = TEST_RETURN_CODE_FAIL;
    }
    else if (confirm_matches_expected_str_arrays(fixture->expected_input_var_names, actual_var_names, var_count)) {
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

    bmi_status = fixture->bmi_model->get_output_var_names(fixture->bmi_model, actual_var_names);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code getting output variable names");
        result = TEST_RETURN_CODE_FAIL;
    }
    else if (confirm_matches_expected_str_arrays(fixture->expected_output_var_names, actual_var_names, var_count)) {
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
    int grid_value, status;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_grid(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], &grid_value);
        if (status == BMI_FAILURE) {
            printf("\nReturned BMI_FAILURE status code getting grid for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_ints(fixture->expected_grid_ids[i], grid_value)) {
            printf("\nGrid value for are different for %s", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_var_itemsize(TestFixture* fixture)
{
    int item_size, expected_size, status;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_itemsize(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], &item_size);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting item size for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // Only SURF_RUNOFF_SCHEME is int, and the rest (output and input) are doubles
        expected_size = strcmp(fixture->expected_output_and_input_var_names[i], "SURF_RUNOFF_SCHEME") == 0 ? sizeof(int) : sizeof(double);

        if (!confirm_matches_expected_ints(expected_size, item_size)) {
            printf("\nSize for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_var_location(TestFixture* fixture)
{
    int status;
    char actual_value[BMI_MAX_VAR_NAME];
    char* expected;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_location(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], &actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting location for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // Only SURF_RUNOFF_SCHEME is "none", and the rest (output and input) are "node"
        if (strcmp(fixture->expected_output_and_input_var_names[i], "SURF_RUNOFF_SCHEME") == 0)
            expected = "none";
        else
            expected = "node";

        if (!confirm_matches_expected_strs(expected, actual_value)) {
            printf("\nLocation for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}


int test_get_var_units(TestFixture* fixture)
{
    int status;
    char actual_value[BMI_MAX_VAR_NAME];

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

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_units(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], &actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting units for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_strs(expected_unit_values[i], actual_value)) {
            printf("\nUnits for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_var_type(TestFixture* fixture)
{
    int status;
    char actual_value[BMI_MAX_TYPE_NAME];
    char* expected_type;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_type(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], &actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting type for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // Only SURF_RUNOFF_SCHEME is int, and the rest (output and input) are doubles
        if (strcmp(fixture->expected_output_and_input_var_names[i], "SURF_RUNOFF_SCHEME") == 0)
            expected_type = "int";
        else
            expected_type = "double";

        if (!confirm_matches_expected_strs(expected_type, actual_value)) {
            printf("\nType for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_var_nbytes(TestFixture* fixture)
{
    int item_nbytes, expected_nbytes, status;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_nbytes(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], &item_nbytes);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting nbytes for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // Only SURF_RUNOFF_SCHEME is int, and the rest (output and input) are doubles
        // And everything should be of single element, so just the size
        if (strcmp(fixture->expected_output_and_input_var_names[i], "SURF_RUNOFF_SCHEME") == 0)
            expected_nbytes = sizeof(int);
        else
            expected_nbytes = sizeof(double);

        if (!confirm_matches_expected_ints(expected_nbytes, item_nbytes)){
            printf("\nnbytes for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
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