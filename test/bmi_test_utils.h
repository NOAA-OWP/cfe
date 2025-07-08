#ifndef BMI_TEST_UTILS
#define BMI_TEST_UTILS

#include "bmi_cfe.h"

/* Define some expected constants */
#define EXPECTED_COMPONENT_NAME "The CFE Model"
#define EXPECTED_INPUT_VAR_COUNT 5
#define EXPECTED_OUTPUT_VAR_COUNT 15
#define EXPECTED_TOTAL_VAR_COUNT 20  // This should equal the sum of the input and output counts

typedef struct TestFixture
{
    Bmi *bmi_model;
    const char* cfg_file;

} TestFixture;

TestFixture* setup(const char* cfg_file);

void teardown(TestFixture* fixture);

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
char** get_all_bmi_variable_names(Bmi* bmi_model, int* output_var_count, int* input_var_count);

/**
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
void setup_expected_grid_ids(int* grid_id_array);

#endif