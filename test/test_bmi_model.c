#include "general_test_utils.h"
#include "bmi_test_utils.h"
#include "bmi_cfe.h"
#include "cfe.h"

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

int test_get_current_time(TestFixture* fixture)
{
    int bmi_status;
    double current_time;

    // For this, we need to be able to initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test getting current time)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Start by checking the initial value, which should be the same as the start time (0.0)
    bmi_status = fixture->bmi_model->get_current_time(fixture->bmi_model, &current_time);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get initial current time");
        return TEST_RETURN_CODE_FAIL;
    }
    if (!confirm_matches_expected_doubles(0.0, current_time)) {
        printf("\nDid not match initial expected module current time");
        return TEST_RETURN_CODE_FAIL;
    }

    // Next, update the model one time step (failing if we couldn't set the variables properly)
    if (!set_arbitrary_input_variables_before_update(fixture, current_time)) {
        printf("\nCouldn't set module input variables before running first update (while testing getting current time");
        return TEST_RETURN_CODE_FAIL;
    }
    bmi_status = fixture->bmi_model->update(fixture->bmi_model);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to advance module with 'update' to test current time");
        return TEST_RETURN_CODE_FAIL;
    }

    // Now make sure the current time value updates as expected
    bmi_status = fixture->bmi_model->get_current_time(fixture->bmi_model, &current_time);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get current time after advancing 1 time step");
        return TEST_RETURN_CODE_FAIL;
    }
    if (!confirm_matches_expected_doubles(3600.0, current_time)) {
        printf("\nDid not match expected module current time after advancing 1 time step");
        return TEST_RETURN_CODE_FAIL;
    }

    return TEST_RETURN_CODE_PASS;
}

int test_get_end_time(TestFixture* fixture)
{

    // The tests use the standard boolean macros from <stdbool.h>, but CFE also defines and uses its own
    // custom TRUE and FALSE macros, in particular in get_end_time.  So, sanity check them to be safe
#ifdef TRUE
    if (TRUE != true) {
        printf("\nConflicting definitions for stdbool.h 'true' (%i) and CFE internal 'TRUE' (%i) macros", true, TRUE);
        printf("\nThis may have implications beyond this test, but is problematic in particular for get_end_time");
        return TEST_RETURN_CODE_FAIL;
    }
#endif
#ifdef FALSE
    if (FALSE != false) {
        printf("\nConflicting definitions for stdbool.h 'false' (%i) and CFE internal 'FALSE' (%i) macros", false, FALSE);
        printf("\nThis may have implications beyond this test, but is problematic in particular for get_end_time");
        return TEST_RETURN_CODE_FAIL;
    }
#endif

    int bmi_status;
    double end_time;

    // For this, we need to be able to initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test getting end time)");
        return TEST_RETURN_CODE_FAIL;
    }

    bmi_status = fixture->bmi_model->get_end_time(fixture->bmi_model, &end_time);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get end time");
        return TEST_RETURN_CODE_FAIL;
    }

    // The testing for this is odd ... for now, fail with an error in this scenario, and note test needs extension
    if (((cfe_state_struct*)(fixture->bmi_model->data))->is_forcing_from_bmi != true) {
        printf("\nTesting of get_end_time not implemented for if module directly reads its own forcings");
        return TEST_RETURN_CODE_FAIL;
    }

    // Assume expected end time is FLT_MAX
    if (!confirm_matches_expected_doubles(FLT_MAX, end_time)) {
        printf("\nDid not match expected module end time");
        return TEST_RETURN_CODE_FAIL;
    }
    return TEST_RETURN_CODE_PASS;
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

int test_get_grid_size(TestFixture* fixture)
{
    int status, actual_size;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_grid_rank(fixture->bmi_model, fixture->expected_grid_ids[i], &actual_size);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting grid size for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // For now, all grids are of size 1
        if (!confirm_matches_expected_ints(1, actual_size)) {
            printf("\nGrid size for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_type(TestFixture* fixture)
{
    int status;

    char actual_type[BMI_MAX_COMPONENT_NAME];

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_grid_type(fixture->bmi_model, fixture->expected_grid_ids[i], actual_type);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting grid type for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // For now, all grids are "scalar"
        if (!confirm_matches_expected_strs("scalar", actual_type)) {
            printf("\nGrid type for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
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

int test_get_start_time(TestFixture* fixture)
{
    int bmi_status;
    double start_time;
    bmi_status = fixture->bmi_model->get_start_time(fixture->bmi_model, &start_time);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get start time");
        return TEST_RETURN_CODE_FAIL;
    }
    if (!confirm_matches_expected_doubles(EXPECTED_MODULE_START_TIME, start_time)) {
        printf("\nDid not match expected module start time");
        return TEST_RETURN_CODE_FAIL;
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_time_step(TestFixture* fixture)
{
    int bmi_status;
    double time_step;
    bmi_status = fixture->bmi_model->get_time_step(fixture->bmi_model, &time_step);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get time step");
        return TEST_RETURN_CODE_FAIL;
    }
    if (!confirm_matches_expected_doubles(EXPECTED_TIME_STEP_SIZE, time_step)) {
        printf("\nDid not match expected module time step");
        return TEST_RETURN_CODE_FAIL;
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_time_units(TestFixture* fixture)
{
    int bmi_status;
    char time_units[BMI_MAX_UNITS_NAME];
    bmi_status = fixture->bmi_model->get_time_units(fixture->bmi_model, time_units);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE status code attempting to get time units");
        return TEST_RETURN_CODE_FAIL;
    }
    // Assume expected time units of seconds (s)
    if (!confirm_matches_expected_strs("s", time_units)) {
        printf("\nDid not match expected module time step");
        return TEST_RETURN_CODE_FAIL;
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_value(TestFixture* fixture)
{
    // For this, we need to be able to initialize first
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_value)");
        return TEST_RETURN_CODE_FAIL;
    }

    char var_type[BMI_MAX_TYPE_NAME];
    double uninit_value, var_value;

    double arbitrary_values[EXPECTED_INPUT_VAR_COUNT];
    get_arbitrary_input_var_values(fixture->current_test_example, 0, arbitrary_values);

    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++) {
        // Have local var for these just for readability
        const char* var_name = fixture->expected_input_var_names[i];
        double* current_arb_val = arbitrary_values + i;

        // Sanity check the test's validity
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, var_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code checking type for '%s' (while testing get_value)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_strs("double", var_type)) {
            printf("\nUnexpected variable type for '%s' (while testing get_value)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Confirm that uninitialized values are not the same thing we will be setting and then getting for the tests
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, &uninit_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting uninitialized value for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        // The odds of this are very low, but we are better off knowing if this happens (though wait until then to
        // do anything about it)
        if (*current_arb_val == uninit_value) {
            printf("\nUninitialized value for '%s' matches test value for get_value and invalidates test", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Then set the arbitrary value and test getting it
        bmi_status = fixture->bmi_model->set_value(fixture->bmi_model, var_name, current_arb_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status setting arbitrary value (for get_value) for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, &var_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting arbitrary value for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(*current_arb_val, var_value)) {
            printf("\nArbitrary value retrieved with get_value was not as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    // Now do some checking of the output variables

    // All doubles except one int, which we will cast for storage purpose
    double output_double_vars[EXPECTED_OUTPUT_VAR_COUNT];
    if (!get_output_var_values(fixture, output_double_vars)) {
        printf("\nFailed to get output var values for testing get_value");
        return TEST_RETURN_CODE_FAIL;
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_value_at_indices(TestFixture* fixture)
{
    // This should basically work the same as the test for get_value, just setting only at indices 0

    // For this, we need to be able to initialize first
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_value_at_indices)");
        return TEST_RETURN_CODE_FAIL;
    }

    char var_type[BMI_MAX_TYPE_NAME];
    double uninit_value, var_value;
    int indices[1] = {0};

    double arbitrary_values[EXPECTED_INPUT_VAR_COUNT];
    get_arbitrary_input_var_values(fixture->current_test_example, 0, arbitrary_values);

    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++) {
        // Have local var for these just for readability
        const char* var_name = fixture->expected_input_var_names[i];
        double* current_arb_val = arbitrary_values + i;

        // Sanity check the test's validity
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, var_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nBMI_FAILURE status checking type for '%s' (while testing set_value_at_indices)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_strs("double", var_type)) {
            printf("\nUnexpected variable type for '%s' (while testing get_value_at_indices)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Confirm that uninitialized values are not the same thing we will be setting and then getting for the tests
        bmi_status = fixture->bmi_model->get_value_at_indices(fixture->bmi_model, var_name, &uninit_value, indices, 1);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting uninitialized value at index 0 for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        // The odds of this are very low, but we are better off knowing if this happens (though wait until then to
        // do anything about it)
        if (*current_arb_val == uninit_value) {
            printf("\nUninitialized value for '%s' matches test value for get_value_at_indices and invalidates test",
                   var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Then set the arbitrary value for an index and test getting it
        bmi_status = fixture->bmi_model->set_value_at_indices(fixture->bmi_model, var_name, indices, 1, current_arb_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE setting arbitrary value (for get_value_at_indices) at index 0 for '%s'",
                   var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        bmi_status = fixture->bmi_model->get_value_at_indices(fixture->bmi_model, var_name, &var_value, indices, 1);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting arbitrary value at index 0 for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(*current_arb_val, var_value)) {
            printf("\nArbitrary value retrieved with get_value_at_indices for index 0 was not as expected for '%s'",
                   var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_value_ptr(TestFixture* fixture)
{
    // For this, we need to be able to initialize first
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test set_value)");
        return TEST_RETURN_CODE_FAIL;
    }

    char var_type[BMI_MAX_TYPE_NAME];
    double uninit_value, var_value, initial_ptr_set_value;
    double arbitrary_values[EXPECTED_INPUT_VAR_COUNT];
    get_arbitrary_input_var_values(fixture->current_test_example, 0, arbitrary_values);

    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++) {
        // Have local var for these just for readability
        const char* var_name = fixture->expected_input_var_names[i];
        double* current_arb_val = arbitrary_values + i;
        double* var_ptr;
        void** var_ptr_ptr = (void**) &var_ptr;

        // Sanity check the test's validity
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, var_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code checking type for '%s' (while testing get_value_ptr)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_strs("double", var_type)) {
            printf("\nUnexpected variable type for '%s' (while testing get_value_ptr)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Get the pointer
        bmi_status = fixture->bmi_model->get_value_ptr(fixture->bmi_model, var_name, var_ptr_ptr);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting pointer for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, &uninit_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting uninitialized value for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        // The odds of this are low, but we are better off knowing if this happens (though wait until then to
        // do anything about it)
        if (0.0 == uninit_value) {
            printf("\nWARN: uninitialized value for '%s' was 0.0; using 0.1 in first set testing pointer", var_name);
            initial_ptr_set_value = 0.1;
        }
        else
            initial_ptr_set_value = 0.0;

        // Now use the pointer to set, setting zero value, and then confirm via get_value things were set right
        *var_ptr = initial_ptr_set_value;
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, &var_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status getting zero value (for get_value_ptr) for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(initial_ptr_set_value, var_value)) {
            printf("\nZero value was not set as expected via pointer for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Finally, set the arbitrary value and make sure it is reflected in the pointer
        bmi_status = fixture->bmi_model->set_value(fixture->bmi_model, var_name, current_arb_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status from set_value (in test for get_value_ptr) for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(*current_arb_val, *var_ptr)) {
            printf("\nArbitrary value retrieved via pointer was not as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    for (int i = 0; i < EXPECTED_OUTPUT_VAR_COUNT; i++) {
        int* int_ptr;
        double* double_ptr;

        void** ptr = (i == 13) ? &int_ptr : &double_ptr;

        const char* var_name = fixture->expected_output_var_names[i];

        // Sanity check the test's validity
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, var_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code checking type for output '%s' (while testing get_value_ptr)",
                   var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        bmi_status = fixture->bmi_model->get_value_ptr(fixture->bmi_model, var_name, ptr);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE getting pointer for output '%s' (while testing get_value_ptr)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        int int_var_val;
        double double_var_val;
        void* var_val = (i == 13) ? &int_var_val : &double_var_val;
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, var_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE getting value for output '%s' (while testing get_value_ptr)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (i == 13 && !confirm_matches_expected_ints(int_var_val, *int_ptr)) {
            printf("\nOutput value retrieved via int pointer was not as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (i != 13 && !confirm_matches_expected_doubles(double_var_val, *double_ptr)) {
            printf("\nOutput value retrieved via double pointer was not as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    return TEST_RETURN_CODE_PASS;
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

int test_set_value(TestFixture* fixture)
{
    // For this, we need to be able to initialize first
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test set_value)");
        return TEST_RETURN_CODE_FAIL;
    }

    char var_type[BMI_MAX_TYPE_NAME];
    double var_value, zero_value = 0.0;

    double arbitrary_values[EXPECTED_INPUT_VAR_COUNT];
    get_arbitrary_input_var_values(fixture->current_test_example, 0, arbitrary_values);

    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++) {
        // Have local var for this just for readability
        const char* var_name = fixture->expected_input_var_names[i];
        // Sanity check the test's validity
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, var_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code checking type for '%s' (while testing set_value)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_strs("double", var_type)) {
            printf("\nUnexpected variable type for '%s' (while testing set_value)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // First, test setting input var to 0
        bmi_status = fixture->bmi_model->set_value(fixture->bmi_model, var_name, &zero_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status setting zero value for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        // And confirm
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, &var_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status confirming zero value was set for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(zero_value, var_value)) {
            printf("\nZero value was not set as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Now, set to the arbitrary value
        double* current_arb_val = arbitrary_values + i;
        bmi_status = fixture->bmi_model->set_value(fixture->bmi_model, var_name, current_arb_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status setting arbitrary value for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        // And confirm
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, &var_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status confirming arbitrary value was set for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(*current_arb_val, var_value)) {
            printf("\nArbitrary value was not set as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_set_value_at_indices(TestFixture* fixture)
{
    // This should basically work the same as the test for set_value, just setting only at indices 0

    // For this, we need to be able to initialize first
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test set_value)");
        return TEST_RETURN_CODE_FAIL;
    }

    char var_type[BMI_MAX_TYPE_NAME];
    double var_value, zero_value = 0.0;
    int indices[1] = {0};

    double arbitrary_values[EXPECTED_INPUT_VAR_COUNT];
    get_arbitrary_input_var_values(fixture->current_test_example, 0, arbitrary_values);

    for (int i = 0; i < EXPECTED_INPUT_VAR_COUNT; i++) {
        // Have local var for this just for readability
        const char* var_name = fixture->expected_input_var_names[i];
        // Sanity check the test's validity
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, var_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nBMI_FAILURE status checking type for '%s' (while testing set_value_at_indices)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_strs("double", var_type)) {
            printf("\nUnexpected variable type for '%s' (while testing set_value_at_indices)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // First, test setting input var to 0
        bmi_status = fixture->bmi_model->set_value_at_indices(fixture->bmi_model, var_name, indices, 1, &zero_value);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status setting zero value at index 0 for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        // And confirm
        bmi_status = fixture->bmi_model->get_value_at_indices(fixture->bmi_model, var_name, &var_value, indices, 1);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status confirming zero value was set at index 0 for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(zero_value, var_value)) {
            printf("\nZero value was not set as expected for at index 0 '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        // Now, set to the arbitrary value
        double* current_arb_val = arbitrary_values + i;
        bmi_status = fixture->bmi_model->set_value_at_indices(fixture->bmi_model, var_name, indices, 1, current_arb_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status setting arbitrary value at index 0 for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        // And confirm
        bmi_status = fixture->bmi_model->get_value_at_indices(fixture->bmi_model, var_name, &var_value, indices, 1);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status confirming arbitrary value at index 0 was set for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(*current_arb_val, var_value)) {
            printf("\nArbitrary value was not set as expected at index 0 for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int main(int argc, const char* argv[])
{
    char* config_file;
    int result = -1;
    unsigned int example_case;

    // Do this sanity check also
    if (EXPECTED_INPUT_VAR_COUNT + EXPECTED_OUTPUT_VAR_COUNT != EXPECTED_TOTAL_VAR_COUNT) {
        printf("\nExpected number of input and output BMI variable constants not consistent with expected total\n");
        printf("\nThose constants are set at %i, %i, and %i respectively",
               EXPECTED_INPUT_VAR_COUNT, EXPECTED_OUTPUT_VAR_COUNT, EXPECTED_TOTAL_VAR_COUNT);
        return TEST_RETURN_CODE_FAIL;
    }

    // Test function is always argv[1]
    // Example test instance is argv[2] (or assumed to be "1" if not given)

    if (argc < 3 || strcmp(argv[2], "1") == 0) {
        config_file = BMI_INIT_CONFIG_EX_1;
        example_case = 1;
    }
    else {
        printf("\nUnexpected test case %s\n", argv[2]);
        return 1;
    }

    TestFixture* fixture = setup(example_case, config_file);

    if (strcmp(argv[1], "test_initialize") == 0)
        result = test_initialize(fixture);
    // TODO: update test
    // TODO: update_until test
    // TODO: finalize test
    else if (strcmp(argv[1], "test_get_component_name") == 0)
        result = test_get_component_name(fixture);
    else if (strcmp(argv[1], "test_get_current_time") == 0)
        result = test_get_current_time(fixture);
    else if (strcmp(argv[1], "test_get_end_time") == 0)
        result = test_get_end_time(fixture);
    else if (strcmp(argv[1], "test_get_grid_rank") == 0)
        result = test_get_grid_rank(fixture);
    else if (strcmp(argv[1], "test_get_grid_size") == 0)
        result = test_get_grid_size(fixture);
    else if (strcmp(argv[1], "test_get_grid_type") == 0)
        result = test_get_grid_type(fixture);
    else if (strcmp(argv[1], "test_get_input_item_count") == 0)
        result = test_get_input_item_count(fixture);
    else if (strcmp(argv[1], "test_get_input_var_names") == 0)
        result = test_get_input_var_names(fixture);
    else if (strcmp(argv[1], "test_get_output_item_count") == 0)
        result = test_get_output_item_count(fixture);
    else if (strcmp(argv[1], "test_get_output_var_names") == 0)
        result = test_get_output_var_names(fixture);
    else if (strcmp(argv[1], "test_get_start_time") == 0)
        result = test_get_start_time(fixture);
    else if (strcmp(argv[1], "test_get_time_step") == 0)
        result = test_get_time_step(fixture);
    else if (strcmp(argv[1], "test_get_time_units") == 0)
        result = test_get_time_units(fixture);
    else if (strcmp(argv[1], "test_get_value") == 0)
        result = test_get_value(fixture);
    else if (strcmp(argv[1], "test_get_value_at_indices") == 0)
        result = test_get_value_at_indices(fixture);
    else if (strcmp(argv[1], "test_get_value_ptr") == 0)
        result = test_get_value_ptr(fixture);
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
    else if (strcmp(argv[1], "test_set_value") == 0)
        result = test_set_value(fixture);
    else if (strcmp(argv[1], "test_set_value_at_indices") == 0)
        result = test_set_value_at_indices(fixture);
    else
        printf("\nUnexpected test function %s\n", argv[1]);

    teardown(fixture);
    return result;
}