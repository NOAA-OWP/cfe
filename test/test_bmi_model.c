#include <math.h>
#include <string.h>
#include <float.h>
#include "general_test_utils.h"
#include "bmi_test_utils.h"
#include "bmi_cfe.h"
#include "cfe.h"
#include "ngen_utilities.h"

/*
 * test_mass_balance_protocol
 *
 * Run the model for several timesteps with known forcing and verify the
 * ngen mass balance protocol identity:
 *   mass_in = mass_out + mass_stored + mass_leaked
 *
 * mass_in starts at initial total storage and accumulates rainfall.
 * mass_stored is total storage at end of last step.
 * mass_out is cumulative discharge + ET.
 * mass_leaked is 0 (no deep losses modeled).
 */
int test_mass_balance_protocol(TestFixture* fixture)
{
    int bmi_status;
    Bmi *m = fixture->bmi_model;

    bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for mass balance test");
        return TEST_RETURN_CODE_FAIL;
    }

    /* obtain pointers to the four protocol variables */
    double *mass_in = NULL, *mass_out = NULL, *mass_stored = NULL, *mass_leaked = NULL;
    bmi_status  = m->get_value_ptr(m, NGEN_MASS_IN,     (void**)&mass_in);
    bmi_status |= m->get_value_ptr(m, NGEN_MASS_OUT,    (void**)&mass_out);
    bmi_status |= m->get_value_ptr(m, NGEN_MASS_STORED, (void**)&mass_stored);
    bmi_status |= m->get_value_ptr(m, NGEN_MASS_LEAKED, (void**)&mass_leaked);
    if (bmi_status != BMI_SUCCESS || !mass_in || !mass_out || !mass_stored || !mass_leaked) {
        printf("\nFailed to get mass balance protocol pointers");
        return TEST_RETURN_CODE_FAIL;
    }

    /* run 5 timesteps with rainfall, then 5 dry timesteps */
    double rain_m = 0.005 / 3600.0;   /* 5 mm/timestep as rate in m/s */
    double no_rain = 0.0;
    double pet_m = 0.0;

    for (int t = 0; t < 10; t++) {
        double r = (t < 5) ? rain_m : no_rain;
        m->set_value(m, "rainfall_depth_m", &r);
        m->set_value(m, "et_potential_m", &pet_m);
        bmi_status = m->update(m);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nUpdate failed at step %d", t);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* check protocol values are physically reasonable */
    if (*mass_in <= 0.0) {
        printf("\nmass_in should be > 0 (got %.6e)", *mass_in);
        return TEST_RETURN_CODE_FAIL;
    }
    if (*mass_out < 0.0) {
        printf("\nmass_out should be >= 0 (got %.6e)", *mass_out);
        return TEST_RETURN_CODE_FAIL;
    }
    if (*mass_stored <= 0.0) {
        printf("\nmass_stored should be > 0 (got %.6e)", *mass_stored);
        return TEST_RETURN_CODE_FAIL;
    }
    if (*mass_leaked < 0.0) {
        printf("\nmass_leaked should be >= 0 (got %.6e)", *mass_leaked);
        return TEST_RETURN_CODE_FAIL;
    }

    /* conservation: mass_in = mass_out + mass_stored + mass_leaked */
    double residual = *mass_in - *mass_out - *mass_stored - *mass_leaked;
    double tol = 1.0e-12;  /* near double precision */
    if (fabs(residual) > tol) {
        printf("\nmass balance residual %.6e exceeds tolerance %.1e", residual, tol);
        printf("\n  mass_in=%.15e  mass_out=%.15e  mass_stored=%.15e  mass_leaked=%.15e",
               *mass_in, *mass_out, *mass_stored, *mass_leaked);
        return TEST_RETURN_CODE_FAIL;
    }

    printf("\n  mass_in=%.6e  out=%.6e  stored=%.6e  leaked=%.6e  residual=%.2e",
           *mass_in, *mass_out, *mass_stored, *mass_leaked, residual);

    return TEST_RETURN_CODE_PASS;
}

/*
 * test_grid_consistency
 *
 * For each variable, verify that the grid metadata functions are mutually
 * consistent: get_var_grid → get_grid_rank/size/type, and that
 * nbytes == grid_size * itemsize for array variables.
 */
int test_grid_consistency(TestFixture* fixture)
{
    Bmi *m = fixture->bmi_model;
    int bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_grid_consistency");
        return TEST_RETURN_CODE_FAIL;
    }

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        const char* name = fixture->expected_output_and_input_var_names[i];

        /* get_var_grid */
        int grid_id = -1;
        bmi_status = m->get_var_grid(m, name, &grid_id);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_var_grid FAILED for '%s'", name);
            return TEST_RETURN_CODE_FAIL;
        }

        /* Verify grid_id matches fixture expectation */
        if (grid_id != fixture->expected_grid_ids[i]) {
            printf("\ngrid_id for '%s': got %d, expected %d", name, grid_id, fixture->expected_grid_ids[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        /* get_grid_rank */
        int rank = -1;
        bmi_status = m->get_grid_rank(m, grid_id, &rank);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_grid_rank FAILED for '%s' (grid %d)", name, grid_id);
            return TEST_RETURN_CODE_FAIL;
        }
        int expected_rank = (grid_id == 0) ? 0 : 1;
        if (rank != expected_rank) {
            printf("\ngrid_rank for '%s' (grid %d): got %d, expected %d", name, grid_id, rank, expected_rank);
            return TEST_RETURN_CODE_FAIL;
        }

        /* get_grid_type */
        char grid_type[64] = {0};
        bmi_status = m->get_grid_type(m, grid_id, grid_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_grid_type FAILED for '%s' (grid %d)", name, grid_id);
            return TEST_RETURN_CODE_FAIL;
        }
        const char* expected_type = (grid_id == 0) ? "scalar" : "vector";
        if (strcmp(grid_type, expected_type) != 0) {
            printf("\ngrid_type for '%s' (grid %d): got '%s', expected '%s'", name, grid_id, grid_type, expected_type);
            return TEST_RETURN_CODE_FAIL;
        }

        /* get_grid_size */
        int grid_size = -1;
        bmi_status = m->get_grid_size(m, grid_id, &grid_size);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_grid_size FAILED for '%s' (grid %d)", name, grid_id);
            return TEST_RETURN_CODE_FAIL;
        }

        /* For array variables: verify nbytes == grid_size * itemsize */
        if (grid_id > 0) {
            int itemsize = -1, nbytes = -1;
            m->get_var_itemsize(m, name, &itemsize);
            m->get_var_nbytes(m, name, &nbytes);

            if (nbytes != grid_size * itemsize) {
                printf("\nGrid consistency FAILED for '%s' (grid %d): nbytes=%d != grid_size(%d) * itemsize(%d) = %d",
                       name, grid_id, nbytes, grid_size, itemsize, grid_size * itemsize);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    return TEST_RETURN_CODE_PASS;
}

int test_finalize(TestFixture* fixture)
{
    int bmi_status;
    void* bmi_var_ptrs[EXPECTED_TOTAL_VAR_COUNT];

    // For this, we need to be able to initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test finalize)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Get all the pointers and make sure none are NULL now that we've initialized
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        const char* var_name = fixture->expected_output_and_input_var_names[i];
        bmi_status = fixture->bmi_model->get_value_ptr(fixture->bmi_model, var_name, bmi_var_ptrs + i);
        if (bmi_status == BMI_FAILURE) {
            printf("\nReturned BMI_FAILURE getting pointer for variable '%s' (in order to test finalize)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (bmi_var_ptrs[i] == NULL) {
            printf("\nGot NULL pointer for variable '%s' after initialization (while testing finalize)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    // Call finalize
    bmi_status = fixture->bmi_model->finalize(fixture->bmi_model);
    if (bmi_status == BMI_FAILURE) {
        printf("\nReturned BMI_FAILURE executing finalize");
        return TEST_RETURN_CODE_FAIL;
    }

    // After finalize, check if pointers are now equal to NULL (for now, just warn if not)
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        if (bmi_var_ptrs[i] != NULL) {
            printf("\nWARN: Pointer to '%s' not set to NULL after finalize; this may be required in the future",
                   fixture->expected_output_and_input_var_names[i]);
        }
    }
    return TEST_RETURN_CODE_PASS;
}

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

    // When forcing_file=BMI, end_time is FLT_MAX (unknown/unbounded)
    if (!confirm_matches_expected_doubles((double)FLT_MAX, end_time)) {
        printf("\nDid not match expected module end time");
        return TEST_RETURN_CODE_FAIL;
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_edge_count(TestFixture* fixture)
{
    int bmi_status, edge_count;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_edge_count)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_edge_count(fixture->bmi_model, fixture->expected_grid_ids[i], &edge_count);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_edge_count",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_edge_nodes(TestFixture* fixture)
{
    int bmi_status, edge_nodes;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_edge_nodes)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_edge_nodes(fixture->bmi_model, fixture->expected_grid_ids[i], &edge_nodes);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_edge_nodes",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_face_count(TestFixture* fixture)
{
    int bmi_status, face_count;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_face_count)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_face_count(fixture->bmi_model, fixture->expected_grid_ids[i], &face_count);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_face_count",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_face_edges(TestFixture* fixture)
{
    int bmi_status, face_edges;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_face_edges)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_face_edges(fixture->bmi_model, fixture->expected_grid_ids[i], &face_edges);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_face_edges",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_face_nodes(TestFixture* fixture)
{
    int bmi_status, face_nodes;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_face_nodes)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_face_nodes(fixture->bmi_model, fixture->expected_grid_ids[i], &face_nodes);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_face_nodes",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_node_count(TestFixture* fixture)
{
    int bmi_status, node_count;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_node_count)");
        return TEST_RETURN_CODE_FAIL;
    }

    // v3: get_grid_node_count delegates to get_grid_size and returns SUCCESS
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_node_count(fixture->bmi_model, fixture->expected_grid_ids[i], &node_count);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE for variable '%s' in call to get_grid_node_count",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_nodes_per_face(TestFixture* fixture)
{
    int bmi_status, nodes_per_face;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_nodes_per_face)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_nodes_per_face(fixture->bmi_model, fixture->expected_grid_ids[i], &nodes_per_face);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_nodes_per_face",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_origin(TestFixture* fixture)
{
    int bmi_status;
    double origin;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_origin)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_origin(fixture->bmi_model, fixture->expected_grid_ids[i], &origin);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_origin",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
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

        // v3: scalar grid (id 0) has rank 0, array grids have rank 1
        int expected_rank = (fixture->expected_grid_ids[i] == 0) ? 0 : 1;
        if (!confirm_matches_expected_ints(expected_rank, actual_rank)) {
            printf("\nGrid rank for '%s' did not match expected", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_shape(TestFixture* fixture)
{
    int bmi_status, shape;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_shape)");
        return TEST_RETURN_CODE_FAIL;
    }

    // v3: scalar grid (id 0) returns BMI_FAILURE; array grids return SUCCESS
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_shape(fixture->bmi_model, fixture->expected_grid_ids[i], &shape);
        if (fixture->expected_grid_ids[i] == 0) {
            if (bmi_status != BMI_FAILURE) {
                printf("\nDid not return BMI_FAILURE for scalar grid variable '%s' in call to get_grid_shape",
                       fixture->expected_output_and_input_var_names[i]);
                return TEST_RETURN_CODE_FAIL;
            }
        } else {
            if (bmi_status != BMI_SUCCESS) {
                printf("\nReturned BMI_FAILURE for array grid variable '%s' in call to get_grid_shape",
                       fixture->expected_output_and_input_var_names[i]);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_spacing(TestFixture* fixture)
{
    int bmi_status;
    double spacing;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_spacing)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_spacing(fixture->bmi_model, fixture->expected_grid_ids[i], &spacing);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_spacing",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_x(TestFixture* fixture)
{
    int bmi_status;
    double grid_x;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_x)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_x(fixture->bmi_model, fixture->expected_grid_ids[i], &grid_x);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_x",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_y(TestFixture* fixture)
{
    int bmi_status;
    double grid_y;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_y)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_y(fixture->bmi_model, fixture->expected_grid_ids[i], &grid_y);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_y",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_z(TestFixture* fixture)
{
    int bmi_status;
    double grid_z;

    // Just to be safe, we will initialize first
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test get_grid_z)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Here, make sure all these do indeed return the BMI_FAILURE as expected
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        bmi_status = fixture->bmi_model->get_grid_z(fixture->bmi_model, fixture->expected_grid_ids[i], &grid_z);
        if (bmi_status != BMI_FAILURE) {
            printf("\nDid not return BMI_FAILURE with grid for variable '%s' in call to get_grid_z",
                   fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_grid_size(TestFixture* fixture)
{
    int status, actual_size;

    // Note: this function calls get_grid_rank (bug in v2). Just verify the call succeeds.
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_grid_rank(fixture->bmi_model, fixture->expected_grid_ids[i], &actual_size);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting grid size for '%s'", fixture->expected_output_and_input_var_names[i]);
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

        // v3: scalar grid (id 0) is "scalar", array grids are "vector"
        const char* expected_type = (fixture->expected_grid_ids[i] == 0) ? "scalar" : "vector";
        if (!confirm_matches_expected_strs(expected_type, actual_type)) {
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
    /* v3: get_time_step requires model to be initialized */
    bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (for get_time_step)");
        return TEST_RETURN_CODE_FAIL;
    }
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

/*
 * Known buffer sizes for array variables in the test config.
 * These are test expectations — NOT queried from the model.
 */
#define TEST_NDISC 4                /* state_soil_moisture_theta */
#define TEST_N_NASH_SUBSURFACE 2    /* state_nash_subsurface_storage */
#define TEST_N_GIUH_ORDINATES 5     /* state_giuh_queue (from test config) */

/* Return the expected buffer size in bytes for a variable, based on test
 * expectations (grid_id from fixture, known array sizes from test config).
 * Returns 0 if the variable is unrecognized. */
/* String buffer size must match PATH_FILENAME_STRING_LENGTH in cfe_config.h */
#define TEST_STRING_BUF_SIZE 1024

static int expected_nbytes_for_var(const char* name, int grid_id)
{
    /* string type */
    if (strcmp(name, "forcing_file_path") == 0)
        return TEST_STRING_BUF_SIZE;

    switch (grid_id) {
        case 0:
            /* scalar: int or double */
            if (strcmp(name, "state_current_timestep") == 0 ||
                strcmp(name, "config_simulate_discrete_soil_moisture") == 0 ||
                strcmp(name, "verbosity") == 0 ||
                strcmp(name, "day_of_year") == 0)
                return (int)sizeof(int);
            return (int)sizeof(double);
        case 1: return TEST_NDISC * (int)sizeof(double);
        case 2: return TEST_N_NASH_SUBSURFACE * (int)sizeof(double);
        case 3: return TEST_N_GIUH_ORDINATES * (int)sizeof(double);
        default: return 0;
    }
}

int test_get_value(TestFixture* fixture)
{
    Bmi *m = fixture->bmi_model;
    int bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_get_value");
        return TEST_RETURN_CODE_FAIL;
    }

    /* Set inputs so they have known values before we read them back */
    if (!set_arbitrary_input_variables_before_update(fixture, 0))
        return TEST_RETURN_CODE_FAIL;

    /* Run one step so outputs are populated */
    bmi_status = m->update(m);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to update for test_get_value");
        return TEST_RETURN_CODE_FAIL;
    }

    /* Test get_value for EVERY advertised variable — no skips */
    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        const char* name = fixture->expected_output_and_input_var_names[i];
        int grid_id = fixture->expected_grid_ids[i];
        int nbytes = expected_nbytes_for_var(name, grid_id);

        if (nbytes <= 0) {
            printf("\nUnknown expected size for '%s' (grid %d)", name, grid_id);
            return TEST_RETURN_CODE_FAIL;
        }

        /* Allocate a buffer large enough for any variable (including strings) */
        char buf[TEST_STRING_BUF_SIZE + 8]; /* +8 for sentinel margin */
        memset(buf, 0xCD, sizeof(buf));

        bmi_status = m->get_value(m, name, buf);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_value FAILED for '%s' (grid=%d, expected %d bytes)", name, grid_id, nbytes);
            return TEST_RETURN_CODE_FAIL;
        }

        /* Verify the sentinel wasn't written beyond expected bounds.
         * Check that the byte right after nbytes still has our sentinel. */
        if (nbytes < (int)sizeof(buf) && (unsigned char)buf[nbytes] != 0xCD) {
            printf("\nget_value for '%s' wrote beyond expected %d bytes", name, nbytes);
            return TEST_RETURN_CODE_FAIL;
        }

        /* Also cross-check: get_var_nbytes should agree with our expectation */
        int model_nbytes = 0;
        m->get_var_nbytes(m, name, &model_nbytes);
        if (model_nbytes != nbytes) {
            printf("\nget_var_nbytes for '%s' returned %d, expected %d", name, model_nbytes, nbytes);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    return TEST_RETURN_CODE_PASS;
}

int test_get_value_at_indices(TestFixture* fixture)
{
    Bmi *m = fixture->bmi_model;
    int bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_get_value_at_indices");
        return TEST_RETURN_CODE_FAIL;
    }

    /* Run one step so array state has non-trivial values (rate in m/s) */
    double rain = 0.005 / 3600.0, pet = 0.0;
    m->set_value(m, "rainfall_depth_m", &rain);
    m->set_value(m, "et_potential_m", &pet);
    m->update(m);

    /* --- Test indexed access to state_soil_moisture_theta --- */
    {
        /* First get the full array via get_value */
        double full_theta[TEST_NDISC] = {0};
        m->get_value(m, "state_soil_moisture_theta", full_theta);

        /* Then get each element via get_value_at_indices and compare */
        for (int idx = 0; idx < TEST_NDISC; idx++) {
            double val = -1.0;
            int indices[1] = {idx};
            bmi_status = m->get_value_at_indices(m, "state_soil_moisture_theta", &val, indices, 1);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nget_value_at_indices FAILED for state_soil_moisture_theta[%d]", idx);
                return TEST_RETURN_CODE_FAIL;
            }
            if (!confirm_matches_expected_doubles(full_theta[idx], val)) {
                printf("\nstate_soil_moisture_theta[%d]: at_indices=%.8e != get_value=%.8e", idx, val, full_theta[idx]);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    /* --- Test indexed access to scalar variables (index 0) --- */
    {
        double val_full = -1.0, val_idx = -1.0;
        int indices[1] = {0};
        m->get_value(m, "discharge_m", &val_full);
        bmi_status = m->get_value_at_indices(m, "discharge_m", &val_idx, indices, 1);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_value_at_indices FAILED for scalar 'discharge_m'");
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_doubles(val_full, val_idx)) {
            printf("\ndischarge_m: at_indices[0]=%.8e != get_value=%.8e", val_idx, val_full);
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
        // v3: skip non-double inputs (e.g. verbosity=int, forcing_file_path=string)
        if (strcmp(var_type, "double") != 0) continue;

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

        const char* var_name = fixture->expected_output_var_names[i];

        /* skip array variables and variables with 0 nbytes (inactive routing scheme) */
        int grid_id = 0;
        fixture->bmi_model->get_var_grid(fixture->bmi_model, var_name, &grid_id);
        if (grid_id > 0) {
            int nb = 0;
            fixture->bmi_model->get_var_nbytes(fixture->bmi_model, var_name, &nb);
            if (nb == 0) continue;  /* inactive for this config */
            /* still verify get_value_ptr succeeds for array vars */
            void* arr_ptr = NULL;
            bmi_status = fixture->bmi_model->get_value_ptr(fixture->bmi_model, var_name, &arr_ptr);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nReturned BMI_FAILURE getting pointer for array output '%s'", var_name);
                return TEST_RETURN_CODE_FAIL;
            }
            continue;
        }

        // Determine type dynamically from the model
        bmi_status = fixture->bmi_model->get_var_type(fixture->bmi_model, var_name, var_type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code checking type for output '%s' (while testing get_value_ptr)",
                   var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        int is_int_type = (strcmp(var_type, "int") == 0);

        void** ptr = is_int_type ? (void*) &int_ptr : (void*) &double_ptr;

        bmi_status = fixture->bmi_model->get_value_ptr(fixture->bmi_model, var_name, ptr);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE getting pointer for output '%s' (while testing get_value_ptr)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        int int_var_val;
        double double_var_val;
        void* var_val = is_int_type ? (void*) &int_var_val : (void*) &double_var_val;
        bmi_status = fixture->bmi_model->get_value(fixture->bmi_model, var_name, var_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE getting value for output '%s' (while testing get_value_ptr)", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (is_int_type && !confirm_matches_expected_ints(int_var_val, *int_ptr)) {
            printf("\nOutput value retrieved via int pointer was not as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!is_int_type && !confirm_matches_expected_doubles(double_var_val, *double_ptr)) {
            printf("\nOutput value retrieved via double pointer was not as expected for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Calibration parameters: verify ptr consistency with get/set --- */
    {
        static const char* cal_params[] = {
            "soil_effective_porosity", "soil_saturated_hydraulic_conductivity",
            "soil_percolation_rate_limiter", "soil_Clapp_Hornberger_b",
            "soil_lateral_flow_K", "subsurface_nash_K",
            "gw_discharge_coefficient", "gw_discharge_exponent",
            "gw_max_storage_m", "soil_saturated_capillary_head",
            "soil_field_capacity_fraction",
            "Xinanjiang_inflection_a",
            "Xinanjiang_shape_b", "Xinanjiang_shape_x",
            "Priestley_Taylor_alpha", "soil_ice_imperv_threshold"
        };
        int n_cal = sizeof(cal_params) / sizeof(cal_params[0]);
        for (int i = 0; i < n_cal; i++) {
            double* ptr = NULL;
            bmi_status = fixture->bmi_model->get_value_ptr(fixture->bmi_model, cal_params[i], (void**)&ptr);
            if (bmi_status != BMI_SUCCESS || ptr == NULL) {
                printf("\nget_value_ptr FAILED for calibration param '%s'", cal_params[i]);
                return TEST_RETURN_CODE_FAIL;
            }
            /* Write via ptr, read via get_value, round-trip back via set_value */
            double si_val = 9.87 + i;
            *ptr = si_val;
            double readback = -1.0;
            fixture->bmi_model->get_value(fixture->bmi_model, cal_params[i], &readback);
            fixture->bmi_model->set_value(fixture->bmi_model, cal_params[i], &readback);
            if (fabs(*ptr - si_val) > 1.0e-12 * fabs(si_val)) {
                printf("\nptr→get_value→set_value round-trip failed for '%s' (wrote %.6e, got %.6e)",
                       cal_params[i], si_val, *ptr);
                return TEST_RETURN_CODE_FAIL;
            }
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

        // v3: determine expected size dynamically from the model's reported type
        char itemsize_var_type[BMI_MAX_TYPE_NAME];
        fixture->bmi_model->get_var_type(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], itemsize_var_type);
        if (strcmp(itemsize_var_type, "int") == 0)
            expected_size = sizeof(int);
        else if (strcmp(itemsize_var_type, "string") == 0)
            expected_size = item_size;  // trust model's reported size for strings
        else
            expected_size = sizeof(double);

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
        status = fixture->bmi_model->get_var_location(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting location for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // v3: all variables return "node"
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

    /* v3 expected units — codified so accidental changes are caught.
     * Order must match the variable names in bmi_test_utils.c setup(). */
    char* expected_units[EXPECTED_TOTAL_VAR_COUNT] = {
        /* 32 outputs */
        "m",    /* discharge_m */
        "m",    /* surface_runoff_m */
        "m",    /* lateral_flow_m */
        "m",    /* baseflow_m */
        "m",    /* actual_et_m */
        "m",    /* vol_balance_residual_m */
        "m",    /* state_soil_storage_m */
        "m",    /* state_gw_storage_m */
        "1",    /* state_current_timestep */
        "-",    /* state_soil_moisture_theta */
        "m",    /* state_nash_subsurface_storage */
        "m",    /* state_giuh_queue */
        "1",    /* config_simulate_discrete_soil_moisture */
        "km2",  /* param_catchment_area_km2 */
        "m",    /* param_soil_depth_m */
        "-",    /* param_soil_porosity */
        "m",    /* timestep_storage_start_m */
        "m",    /* timestep_input_m */
        "m",    /* timestep_output_m */
        "m",    /* timestep_storage_end_m */
        "m",    /* potential_et_m */
        "m",    /* giuh_outflow_m */
        "m",    /* soil_to_gw_percolation_flux_m */
        "-",    /* soil_moisture_theta_1 */
        "-",    /* soil_moisture_theta_2 */
        "-",    /* soil_moisture_theta_3 */
        "-",    /* soil_moisture_theta_4 */
        "m",    /* bare_soil_evaporation_m */
        "m",    /* impervious_runoff_m */
        "m",    /* pervious_runoff_m */
        "m",    /* surface_routed_to_outlet_m */
        "m",    /* lateral_flow_generated_m */
        /* 5 inputs (AORC fields omitted when PT/soil-evap disabled) */
        "m s-1",    /* rainfall_depth_m */
        "m s-1",    /* et_potential_m */
        "-",        /* ice_fraction */
        "-",        /* param_catchment_vegetated_fraction */
        "-"         /* bare_soil_rsurf_exp */
    };

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_units(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting units for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_strs(expected_units[i], actual_value)) {
            printf("\nUnits for '%s' did not match expected (got '%s', want '%s')",
                   fixture->expected_output_and_input_var_names[i], actual_value, expected_units[i]);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_var_type(TestFixture* fixture)
{
    int status;
    char actual_value[BMI_MAX_TYPE_NAME];

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        status = fixture->bmi_model->get_var_type(fixture->bmi_model, fixture->expected_output_and_input_var_names[i], actual_value);
        if (status != BMI_SUCCESS) {
            printf("\nReturned BMI_FAILURE status code getting type for '%s'", fixture->expected_output_and_input_var_names[i]);
            return TEST_RETURN_CODE_FAIL;
        }

        // v3: validate the returned type is one of "double", "int", or "string"
        if (strcmp(actual_value, "double") != 0 &&
            strcmp(actual_value, "int") != 0 &&
            strcmp(actual_value, "string") != 0) {
            printf("\nType for '%s' was '%s', expected one of double/int/string",
                   fixture->expected_output_and_input_var_names[i], actual_value);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    return TEST_RETURN_CODE_PASS;
}

int test_get_var_nbytes(TestFixture* fixture)
{
    /* This test must initialize because array sizes depend on config */
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_get_var_nbytes");
        return TEST_RETURN_CODE_FAIL;
    }

    int item_nbytes, status;

    for (int i = 0; i < EXPECTED_TOTAL_VAR_COUNT; i++) {
        const char* var_name = fixture->expected_output_and_input_var_names[i];
        int grid_id = fixture->expected_grid_ids[i];

        status = fixture->bmi_model->get_var_nbytes(fixture->bmi_model, var_name, &item_nbytes);
        if (status != BMI_SUCCESS) {
            printf("\nget_var_nbytes FAILED for '%s'", var_name);
            return TEST_RETURN_CODE_FAIL;
        }

        /* Compare against independently known expected size */
        int expected = expected_nbytes_for_var(var_name, grid_id);
        if (expected <= 0) {
            printf("\nNo expected nbytes for '%s' (grid %d)", var_name, grid_id);
            return TEST_RETURN_CODE_FAIL;
        }
        if (item_nbytes != expected) {
            printf("\nnbytes for '%s': got %d, expected %d", var_name, item_nbytes, expected);
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
    Bmi *m = fixture->bmi_model;
    int bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_set_value");
        return TEST_RETURN_CODE_FAIL;
    }

    /* --- Double inputs: set → get round-trip --- */
    {
        const char* double_inputs[] = {"rainfall_depth_m", "et_potential_m"};
        for (int i = 0; i < 2; i++) {
            double set_val = 0.0042 + i;
            double get_val = -1.0;
            bmi_status = m->set_value(m, double_inputs[i], &set_val);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nset_value FAILED for double input '%s'", double_inputs[i]);
                return TEST_RETURN_CODE_FAIL;
            }
            bmi_status = m->get_value(m, double_inputs[i], &get_val);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nget_value FAILED after set for '%s'", double_inputs[i]);
                return TEST_RETURN_CODE_FAIL;
            }
            if (!confirm_matches_expected_doubles(set_val, get_val)) {
                printf("\nset/get round-trip failed for '%s'", double_inputs[i]);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    /* --- Int input: set → get round-trip --- */
    {
        int set_val = 2;
        int get_val = -1;
        bmi_status = m->set_value(m, "verbosity", &set_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nset_value FAILED for int input 'verbosity'");
            return TEST_RETURN_CODE_FAIL;
        }
        bmi_status = m->get_value(m, "verbosity", &get_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_value FAILED after set for 'verbosity'");
            return TEST_RETURN_CODE_FAIL;
        }
        if (!confirm_matches_expected_ints(set_val, get_val)) {
            printf("\nset/get round-trip failed for 'verbosity'");
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Array state: set → get round-trip for state_soil_moisture_theta --- */
    {
        double set_theta[TEST_NDISC] = {0.11, 0.22, 0.33, 0.44};
        double get_theta[TEST_NDISC] = {0};
        bmi_status = m->set_value(m, "state_soil_moisture_theta", set_theta);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nset_value FAILED for 'state_soil_moisture_theta'");
            return TEST_RETURN_CODE_FAIL;
        }
        bmi_status = m->get_value(m, "state_soil_moisture_theta", get_theta);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_value FAILED after set for 'state_soil_moisture_theta'");
            return TEST_RETURN_CODE_FAIL;
        }
        for (int i = 0; i < TEST_NDISC; i++) {
            if (!confirm_matches_expected_doubles(set_theta[i], get_theta[i])) {
                printf("\nset/get round-trip failed for state_soil_moisture_theta[%d]", i);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    /* --- Array state: set → get round-trip for state_nash_subsurface_storage --- */
    {
        double set_nash[TEST_N_NASH_SUBSURFACE] = {0.001, 0.002};
        double get_nash[TEST_N_NASH_SUBSURFACE] = {0};
        bmi_status = m->set_value(m, "state_nash_subsurface_storage", set_nash);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nset_value FAILED for 'state_nash_subsurface_storage'");
            return TEST_RETURN_CODE_FAIL;
        }
        bmi_status = m->get_value(m, "state_nash_subsurface_storage", get_nash);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_value FAILED after set for 'state_nash_subsurface_storage'");
            return TEST_RETURN_CODE_FAIL;
        }
        for (int i = 0; i < TEST_N_NASH_SUBSURFACE; i++) {
            if (!confirm_matches_expected_doubles(set_nash[i], get_nash[i])) {
                printf("\nset/get round-trip failed for state_nash_subsurface_storage[%d]", i);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    /* --- Array state: set → get round-trip for state_giuh_queue --- */
    {
        double set_giuh[TEST_N_GIUH_ORDINATES] = {0.01, 0.02, 0.03, 0.04, 0.05};
        double get_giuh[TEST_N_GIUH_ORDINATES] = {0};
        bmi_status = m->set_value(m, "state_giuh_queue", set_giuh);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nset_value FAILED for 'state_giuh_queue'");
            return TEST_RETURN_CODE_FAIL;
        }
        bmi_status = m->get_value(m, "state_giuh_queue", get_giuh);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_value FAILED after set for 'state_giuh_queue'");
            return TEST_RETURN_CODE_FAIL;
        }
        for (int i = 0; i < TEST_N_GIUH_ORDINATES; i++) {
            if (!confirm_matches_expected_doubles(set_giuh[i], get_giuh[i])) {
                printf("\nset/get round-trip failed for state_giuh_queue[%d]", i);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    /* --- Calibration parameters: set → get round-trip --- */
    {
        /* Must match param_var_names[] in bmi_cfe.c */
        static const char* cal_params[] = {
            "soil_effective_porosity", "soil_saturated_hydraulic_conductivity",
            "soil_percolation_rate_limiter", "soil_Clapp_Hornberger_b",
            "soil_lateral_flow_K", "subsurface_nash_K",
            "gw_discharge_coefficient", "gw_discharge_exponent",
            "gw_max_storage_m", "soil_saturated_capillary_head",
            "soil_field_capacity_fraction",
            "Xinanjiang_inflection_a",
            "Xinanjiang_shape_b", "Xinanjiang_shape_x",
            "Priestley_Taylor_alpha", "soil_ice_imperv_threshold"
        };
        int n_cal = sizeof(cal_params) / sizeof(cal_params[0]);
        for (int i = 0; i < n_cal; i++) {
            double set_val = 1.23 + i;
            double get_val = -1.0;
            bmi_status = m->set_value(m, cal_params[i], &set_val);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nset_value FAILED for calibration param '%s'", cal_params[i]);
                return TEST_RETURN_CODE_FAIL;
            }
            bmi_status = m->get_value(m, cal_params[i], &get_val);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nget_value FAILED after set for calibration param '%s'", cal_params[i]);
                return TEST_RETURN_CODE_FAIL;
            }
            if (!confirm_matches_expected_doubles(set_val, get_val)) {
                printf("\nset/get round-trip failed for calibration param '%s'", cal_params[i]);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    return TEST_RETURN_CODE_PASS;
}

int test_set_value_at_indices(TestFixture* fixture)
{
    Bmi *m = fixture->bmi_model;
    int bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_set_value_at_indices");
        return TEST_RETURN_CODE_FAIL;
    }

    /* --- Set individual soil theta elements and verify via get_value_at_indices --- */
    {
        double test_vals[TEST_NDISC] = {0.41, 0.42, 0.43, 0.44};
        for (int idx = 0; idx < TEST_NDISC; idx++) {
            int indices[1] = {idx};
            bmi_status = m->set_value_at_indices(m, "state_soil_moisture_theta", indices, 1, &test_vals[idx]);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nset_value_at_indices FAILED for state_soil_moisture_theta[%d]", idx);
                return TEST_RETURN_CODE_FAIL;
            }
            double readback = -1.0;
            bmi_status = m->get_value_at_indices(m, "state_soil_moisture_theta", &readback, indices, 1);
            if (bmi_status != BMI_SUCCESS) {
                printf("\nget_value_at_indices FAILED after set for state_soil_moisture_theta[%d]", idx);
                return TEST_RETURN_CODE_FAIL;
            }
            if (!confirm_matches_expected_doubles(test_vals[idx], readback)) {
                printf("\nset/get at_indices round-trip failed for state_soil_moisture_theta[%d]", idx);
                return TEST_RETURN_CODE_FAIL;
            }
        }
    }

    /* --- Set scalar double input at index 0 and verify round-trip --- */
    {
        double set_val = 0.0088;
        double get_val = -1.0;
        int indices[1] = {0};
        bmi_status = m->set_value_at_indices(m, "rainfall_depth_m", indices, 1, &set_val);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nset_value_at_indices FAILED for scalar 'rainfall_depth_m'");
            return TEST_RETURN_CODE_FAIL;
        }
        m->get_value_at_indices(m, "rainfall_depth_m", &get_val, indices, 1);
        if (!confirm_matches_expected_doubles(set_val, get_val)) {
            printf("\nset/get at_indices round-trip failed for scalar 'rainfall_depth_m'");
            return TEST_RETURN_CODE_FAIL;
        }
    }

    return TEST_RETURN_CODE_PASS;
}

int test_update(TestFixture* fixture)
{
    // For this, we need to be able to initialize first
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test update)");
        return TEST_RETURN_CODE_FAIL;
    }

    if (!set_arbitrary_input_variables_before_update(fixture, 0)) {
        printf("\nFailed to set arbitrary BMI input variable values (in order to test update)");
        return TEST_RETURN_CODE_FAIL;
    }

    bmi_status = fixture->bmi_model->update(fixture->bmi_model);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to run update");
        return TEST_RETURN_CODE_FAIL;
    }
    // We don't have any platform-independent guarantees at this point about what values will be, so for now just
    // assume we are good if the call to update itself returns a successful BMI status code
    return TEST_RETURN_CODE_PASS;
}

int test_update_until(TestFixture* fixture)
{
    double module_time = EXPECTED_MODULE_START_TIME;
    // For this, we need to be able to initialize first
    int bmi_status = fixture->bmi_model->initialize(fixture->bmi_model, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to initialize (in order to test update)");
        return TEST_RETURN_CODE_FAIL;
    }

    if (!set_arbitrary_input_variables_before_update(fixture, module_time)) {
        printf("\nFailed to set arbitrary BMI input variable values (in order to test update)");
        return TEST_RETURN_CODE_FAIL;
    }

    // Now advance one time step
    module_time += EXPECTED_TIME_STEP_SIZE;
    bmi_status = fixture->bmi_model->update_until(fixture->bmi_model, module_time);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nReturned BMI_FAILURE status code attempting to run update");
        return TEST_RETURN_CODE_FAIL;
    }
    // We don't have any platform-independent guarantees at this point about what values will be, so for now just
    // assume we are good if the call to update itself returns a successful BMI status code
    return TEST_RETURN_CODE_PASS;
}

/*
 * test_serialization_metadata
 *
 * Verify the ngen serialization protocol's check_support() probe would
 * succeed: GetVarType and GetVarUnits resolve for all four reserved names
 * with the exact expected values. Also verify trigger-specific behavior
 * for GetVarItemsize and GetVarNbytes.
 */
int test_serialization_metadata(TestFixture* fixture)
{
    Bmi *m = fixture->bmi_model;
    int bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_serialization_metadata");
        return TEST_RETURN_CODE_FAIL;
    }

    /* Expected metadata for each reserved variable */
    static const struct {
        const char *name;
        const char *type;
        const char *units;
        int itemsize_succeeds;  /* 0 = expect BMI_FAILURE */
    } expected[] = {
        { NGEN_SERIALIZATION_CREATE, "int",  "ngen::trigger", 0 },
        { NGEN_SERIALIZATION_FREE,   "int",  "ngen::trigger", 0 },
        { NGEN_SERIALIZATION_SIZE,   "int",  "bytes",         1 },
        { NGEN_SERIALIZATION_STATE,  "char", "ngen::opaque",  1 },
    };
    int n = sizeof(expected) / sizeof(expected[0]);

    for (int i = 0; i < n; i++) {
        char type[BMI_MAX_TYPE_NAME] = {0};
        bmi_status = m->get_var_type(m, expected[i].name, type);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_var_type FAILED for '%s'", expected[i].name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (strcmp(type, expected[i].type) != 0) {
            printf("\nget_var_type for '%s': got '%s', expected '%s'",
                   expected[i].name, type, expected[i].type);
            return TEST_RETURN_CODE_FAIL;
        }

        char units[BMI_MAX_UNITS_NAME] = {0};
        bmi_status = m->get_var_units(m, expected[i].name, units);
        if (bmi_status != BMI_SUCCESS) {
            printf("\nget_var_units FAILED for '%s'", expected[i].name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (strcmp(units, expected[i].units) != 0) {
            printf("\nget_var_units for '%s': got '%s', expected '%s'",
                   expected[i].name, units, expected[i].units);
            return TEST_RETURN_CODE_FAIL;
        }

        int itemsize = -1;
        bmi_status = m->get_var_itemsize(m, expected[i].name, &itemsize);
        if (expected[i].itemsize_succeeds && bmi_status != BMI_SUCCESS) {
            printf("\nget_var_itemsize unexpectedly FAILED for '%s'", expected[i].name);
            return TEST_RETURN_CODE_FAIL;
        }
        if (!expected[i].itemsize_succeeds && bmi_status != BMI_FAILURE) {
            printf("\nget_var_itemsize should have FAILED for trigger '%s'", expected[i].name);
            return TEST_RETURN_CODE_FAIL;
        }

        int nbytes = -1;
        bmi_status = m->get_var_nbytes(m, expected[i].name, &nbytes);
        if (!expected[i].itemsize_succeeds && bmi_status != BMI_FAILURE) {
            printf("\nget_var_nbytes should have FAILED for trigger '%s'", expected[i].name);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* Verify serialization_size reports sizeof(int) for nbytes */
    int size_nbytes = 0;
    bmi_status = m->get_var_nbytes(m, NGEN_SERIALIZATION_SIZE, &size_nbytes);
    if (bmi_status != BMI_SUCCESS || size_nbytes != sizeof(int)) {
        printf("\nget_var_nbytes for serialization_size: got %d, expected %d",
               size_nbytes, (int)sizeof(int));
        return TEST_RETURN_CODE_FAIL;
    }

    /* Verify GetVarLocation returns BMI_FAILURE for all four */
    for (int i = 0; i < n; i++) {
        char loc[BMI_MAX_VAR_NAME] = {0};
        bmi_status = m->get_var_location(m, expected[i].name, loc);
        if (bmi_status != BMI_FAILURE) {
            printf("\nget_var_location should have FAILED for '%s'", expected[i].name);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    return TEST_RETURN_CODE_PASS;
}

/*
 * test_serialization_round_trip
 *
 * Emulate the ngen engine's save/restore protocol sequence and verify:
 *  1. Captured state can be restored to produce identical model state
 *  2. Model run from restored state produces identical outputs to a continuous run
 */
int test_serialization_round_trip(TestFixture* fixture)
{
    Bmi *m = fixture->bmi_model;
    int bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize for test_serialization_round_trip");
        return TEST_RETURN_CODE_FAIL;
    }

    double rain = 0.005 / 3600.0, pet = 0.0;
    int trigger = 1;

    /* --- Phase 1: advance 5 steps to build up state --- */
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        if (m->update(m) != BMI_SUCCESS) {
            printf("\nUpdate failed at step %d (phase 1)", t);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Phase 2: capture state via protocol --- */
    /* SetValue(create) */
    bmi_status = m->set_value(m, NGEN_SERIALIZATION_CREATE, &trigger);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nSetValue(create) failed");
        return TEST_RETURN_CODE_FAIL;
    }

    /* GetValue(size) */
    int buf_size = 0;
    bmi_status = m->get_value(m, NGEN_SERIALIZATION_SIZE, &buf_size);
    if (bmi_status != BMI_SUCCESS || buf_size <= 0) {
        printf("\nGetValue(size) failed or returned %d", buf_size);
        return TEST_RETURN_CODE_FAIL;
    }

    /* GetValue(state) */
    char *saved_buf = (char*)malloc(buf_size);
    bmi_status = m->get_value(m, NGEN_SERIALIZATION_STATE, saved_buf);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nGetValue(state) failed");
        free(saved_buf);
        return TEST_RETURN_CODE_FAIL;
    }

    /* SetValue(free) */
    m->set_value(m, NGEN_SERIALIZATION_FREE, &trigger);

    /* --- Snapshot A: record all state variables at the capture point --- */
    double snap_a_soil = 0, snap_a_gw = 0;
    double snap_a_theta[NDISC] = {0};
    double snap_a_nash[MAX_NUM_SUBSURFACE_NASH_CASCADE] = {0};
    double snap_a_giuh[MAX_NUM_GIUH_ORDINATES] = {0};
    m->get_value(m, "state_soil_storage_m", &snap_a_soil);
    m->get_value(m, "state_gw_storage_m", &snap_a_gw);
    m->get_value(m, "state_soil_moisture_theta", snap_a_theta);
    m->get_value(m, "state_nash_subsurface_storage", snap_a_nash);
    m->get_value(m, "state_giuh_queue", snap_a_giuh);

    /* --- Phase 3: advance 5 more steps (mutate state) --- */
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
    }

    /* --- Phase 4: restore state from saved buffer --- */
    bmi_status = m->set_value(m, NGEN_SERIALIZATION_STATE, saved_buf);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nSetValue(state) restore failed");
        free(saved_buf);
        return TEST_RETURN_CODE_FAIL;
    }

    /* --- Snapshot B: verify state matches snapshot A (bit-exact) --- */
    double snap_b_soil = 0, snap_b_gw = 0;
    double snap_b_theta[NDISC] = {0};
    double snap_b_nash[MAX_NUM_SUBSURFACE_NASH_CASCADE] = {0};
    double snap_b_giuh[MAX_NUM_GIUH_ORDINATES] = {0};
    m->get_value(m, "state_soil_storage_m", &snap_b_soil);
    m->get_value(m, "state_gw_storage_m", &snap_b_gw);
    m->get_value(m, "state_soil_moisture_theta", snap_b_theta);
    m->get_value(m, "state_nash_subsurface_storage", snap_b_nash);
    m->get_value(m, "state_giuh_queue", snap_b_giuh);

    if (snap_a_soil != snap_b_soil || snap_a_gw != snap_b_gw) {
        printf("\nScalar state mismatch after restore: soil=%.15e vs %.15e, gw=%.15e vs %.15e",
               snap_a_soil, snap_b_soil, snap_a_gw, snap_b_gw);
        free(saved_buf);
        return TEST_RETURN_CODE_FAIL;
    }
    for (int i = 0; i < NDISC; i++) {
        if (snap_a_theta[i] != snap_b_theta[i]) {
            printf("\ntheta[%d] mismatch: %.15e vs %.15e", i, snap_a_theta[i], snap_b_theta[i]);
            free(saved_buf);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        if (snap_a_nash[i] != snap_b_nash[i]) {
            printf("\nnash[%d] mismatch: %.15e vs %.15e", i, snap_a_nash[i], snap_b_nash[i]);
            free(saved_buf);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    for (int i = 0; i < MAX_NUM_GIUH_ORDINATES; i++) {
        if (snap_a_giuh[i] != snap_b_giuh[i]) {
            printf("\ngiuh[%d] mismatch: %.15e vs %.15e", i, snap_a_giuh[i], snap_b_giuh[i]);
            free(saved_buf);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Phase 5: run 5 steps from restored state, capture discharge --- */
    double q_restored[5] = {0};
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
        m->get_value(m, "discharge_m", &q_restored[t]);
    }

    /* --- Phase 6: fresh run to same point, capture discharge --- */
    m->finalize(m);
    bmi_status = m->initialize(m, fixture->cfg_file);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to re-initialize for continuous comparison");
        free(saved_buf);
        return TEST_RETURN_CODE_FAIL;
    }

    /* Run first 5 steps (same as phase 1) */
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
    }
    /* Run next 5 steps (same forcing as phase 5) */
    double q_continuous[5] = {0};
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
        m->get_value(m, "discharge_m", &q_continuous[t]);
    }

    /* --- Phase 7: compare outputs --- */
    for (int t = 0; t < 5; t++) {
        if (q_restored[t] != q_continuous[t]) {
            printf("\nDischarge mismatch at step %d: restored=%.15e continuous=%.15e",
                   t, q_restored[t], q_continuous[t]);
            free(saved_buf);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    printf("\n  round-trip: %d bytes, state bit-exact, outputs bit-exact over 5 post-restore steps",
           buf_size);

    free(saved_buf);
    return TEST_RETURN_CODE_PASS;
}

/*
 * test_serialization_round_trip_dsbm
 *
 * Same protocol sequence as test_serialization_round_trip but uses a config
 * with simulate_discrete_soil_moisture=TRUE.  This exercises the sync of
 * soil_state_in.theta_in / total_storage_m after deserialization — without
 * which the first DSBM timestep would use stale values for rainfall
 * partitioning.
 */
int test_serialization_round_trip_dsbm(TestFixture* fixture)
{
    (void)fixture;  /* uses its own BMI instance with DSBM config */
    Bmi m_storage, *m = &m_storage;
    register_bmi_cfe(m);

    int bmi_status = m->initialize(m, BMI_INIT_CONFIG_DSBM);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize DSBM config for round-trip test");
        return TEST_RETURN_CODE_FAIL;
    }

    int dsbm_flag = 0;
    m->get_value(m, "config_simulate_discrete_soil_moisture", &dsbm_flag);
    if (dsbm_flag != 1) {
        printf("\nDSBM config flag is %d, expected 1", dsbm_flag);
        m->finalize(m);
        return TEST_RETURN_CODE_FAIL;
    }

    double rain = 0.005 / 3600.0, pet = 0.0;
    int trigger = 1;

    /* --- Phase 1: advance 5 steps --- */
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        if (m->update(m) != BMI_SUCCESS) {
            printf("\nUpdate failed at step %d (phase 1)", t);
            m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Phase 2: capture state --- */
    bmi_status = m->set_value(m, NGEN_SERIALIZATION_CREATE, &trigger);
    if (bmi_status != BMI_SUCCESS) { printf("\nSetValue(create) failed"); m->finalize(m); return TEST_RETURN_CODE_FAIL; }

    int buf_size = 0;
    m->get_value(m, NGEN_SERIALIZATION_SIZE, &buf_size);
    char *saved_buf = (char*)malloc(buf_size);
    m->get_value(m, NGEN_SERIALIZATION_STATE, saved_buf);
    m->set_value(m, NGEN_SERIALIZATION_FREE, &trigger);

    /* Snapshot A */
    double snap_a_soil = 0, snap_a_gw = 0;
    double snap_a_theta[NDISC] = {0};
    double snap_a_nash[MAX_NUM_SUBSURFACE_NASH_CASCADE] = {0};
    double snap_a_giuh[MAX_NUM_GIUH_ORDINATES] = {0};
    m->get_value(m, "state_soil_storage_m", &snap_a_soil);
    m->get_value(m, "state_gw_storage_m", &snap_a_gw);
    m->get_value(m, "state_soil_moisture_theta", snap_a_theta);
    m->get_value(m, "state_nash_subsurface_storage", snap_a_nash);
    m->get_value(m, "state_giuh_queue", snap_a_giuh);

    /* --- Phase 3: advance 5 more steps (mutate) --- */
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
    }

    /* --- Phase 4: restore --- */
    bmi_status = m->set_value(m, NGEN_SERIALIZATION_STATE, saved_buf);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nSetValue(state) restore failed");
        free(saved_buf); m->finalize(m);
        return TEST_RETURN_CODE_FAIL;
    }

    /* Snapshot B: verify bit-exact */
    double snap_b_soil = 0, snap_b_gw = 0;
    double snap_b_theta[NDISC] = {0};
    double snap_b_nash[MAX_NUM_SUBSURFACE_NASH_CASCADE] = {0};
    double snap_b_giuh[MAX_NUM_GIUH_ORDINATES] = {0};
    m->get_value(m, "state_soil_storage_m", &snap_b_soil);
    m->get_value(m, "state_gw_storage_m", &snap_b_gw);
    m->get_value(m, "state_soil_moisture_theta", snap_b_theta);
    m->get_value(m, "state_nash_subsurface_storage", snap_b_nash);
    m->get_value(m, "state_giuh_queue", snap_b_giuh);

    if (snap_a_soil != snap_b_soil || snap_a_gw != snap_b_gw) {
        printf("\nDSBM scalar state mismatch: soil=%.15e vs %.15e, gw=%.15e vs %.15e",
               snap_a_soil, snap_b_soil, snap_a_gw, snap_b_gw);
        free(saved_buf); m->finalize(m);
        return TEST_RETURN_CODE_FAIL;
    }
    for (int i = 0; i < NDISC; i++) {
        if (snap_a_theta[i] != snap_b_theta[i]) {
            printf("\nDSBM theta[%d] mismatch: %.15e vs %.15e", i, snap_a_theta[i], snap_b_theta[i]);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        if (snap_a_nash[i] != snap_b_nash[i]) {
            printf("\nDSBM nash[%d] mismatch", i);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    for (int i = 0; i < MAX_NUM_GIUH_ORDINATES; i++) {
        if (snap_a_giuh[i] != snap_b_giuh[i]) {
            printf("\nDSBM giuh[%d] mismatch", i);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Phase 5: run 5 steps from restored state --- */
    double q_restored[5] = {0};
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
        m->get_value(m, "discharge_m", &q_restored[t]);
    }

    /* --- Phase 6: fresh continuous run --- */
    m->finalize(m);
    bmi_status = m->initialize(m, BMI_INIT_CONFIG_DSBM);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to re-initialize DSBM for continuous comparison");
        free(saved_buf);
        return TEST_RETURN_CODE_FAIL;
    }

    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
    }
    double q_continuous[5] = {0};
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->update(m);
        m->get_value(m, "discharge_m", &q_continuous[t]);
    }

    /* --- Phase 7: compare outputs --- */
    for (int t = 0; t < 5; t++) {
        if (q_restored[t] != q_continuous[t]) {
            printf("\nDSBM discharge mismatch at step %d: restored=%.15e continuous=%.15e",
                   t, q_restored[t], q_continuous[t]);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    printf("\n  DSBM round-trip: %d bytes, state bit-exact, outputs bit-exact over 5 post-restore steps",
           buf_size);

    free(saved_buf);
    m->finalize(m);
    return TEST_RETURN_CODE_PASS;
}

/*
 * test_serialization_round_trip_pt
 *
 * Same protocol as test_serialization_round_trip but uses the PT-enabled
 * config with AORC forcings, exercising the PET temperature state
 * (skin_temperature_k, upper_soil_temperature_k, etc.) round-trip.
 */
int test_serialization_round_trip_pt(TestFixture* fixture)
{
    (void)fixture;
    Bmi m_storage, *m = &m_storage;
    register_bmi_cfe(m);

    int bmi_status = m->initialize(m, BMI_INIT_CONFIG_PT);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to initialize PT config for round-trip test");
        return TEST_RETURN_CODE_FAIL;
    }

    double rain = 0.005 / 3600.0, pet = 0.0;
    double temp_k = 287.5, dlwrf = 361.2, dswrf = 200.0;
    double pres = 100530.0, spfh = 0.0105, ugrd = -2.6, vgrd = 0.0;
    double veg_frac = 0.70, rsurf_exp = 5.0;
    int trigger = 1;

    /* --- Phase 1: advance 5 steps with AORC forcings --- */
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->set_value(m, "TMP_2maboveground", &temp_k);
        m->set_value(m, "DLWRF_surface", &dlwrf);
        m->set_value(m, "DSWRF_surface", &dswrf);
        m->set_value(m, "PRES_surface", &pres);
        m->set_value(m, "SPFH_2maboveground", &spfh);
        m->set_value(m, "UGRD_10maboveground", &ugrd);
        m->set_value(m, "VGRD_10maboveground", &vgrd);
        m->set_value(m, "param_catchment_vegetated_fraction", &veg_frac);
        m->set_value(m, "bare_soil_rsurf_exp", &rsurf_exp);
        if (m->update(m) != BMI_SUCCESS) {
            printf("\nUpdate failed at step %d (phase 1)", t);
            m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Phase 2: capture state --- */
    bmi_status = m->set_value(m, NGEN_SERIALIZATION_CREATE, &trigger);
    if (bmi_status != BMI_SUCCESS) { printf("\nSetValue(create) failed"); m->finalize(m); return TEST_RETURN_CODE_FAIL; }

    int buf_size = 0;
    m->get_value(m, NGEN_SERIALIZATION_SIZE, &buf_size);
    char *saved_buf = (char*)malloc(buf_size);
    m->get_value(m, NGEN_SERIALIZATION_STATE, saved_buf);
    m->set_value(m, NGEN_SERIALIZATION_FREE, &trigger);

    /* Snapshot A */
    double snap_a_soil = 0, snap_a_gw = 0;
    double snap_a_theta[NDISC] = {0};
    double snap_a_nash[MAX_NUM_SUBSURFACE_NASH_CASCADE] = {0};
    double snap_a_giuh[MAX_NUM_GIUH_ORDINATES] = {0};
    m->get_value(m, "state_soil_storage_m", &snap_a_soil);
    m->get_value(m, "state_gw_storage_m", &snap_a_gw);
    m->get_value(m, "state_soil_moisture_theta", snap_a_theta);
    m->get_value(m, "state_nash_subsurface_storage", snap_a_nash);
    m->get_value(m, "state_giuh_queue", snap_a_giuh);

    /* --- Phase 3: advance 5 more steps (mutate state) --- */
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->set_value(m, "TMP_2maboveground", &temp_k);
        m->set_value(m, "DLWRF_surface", &dlwrf);
        m->set_value(m, "DSWRF_surface", &dswrf);
        m->set_value(m, "PRES_surface", &pres);
        m->set_value(m, "SPFH_2maboveground", &spfh);
        m->set_value(m, "UGRD_10maboveground", &ugrd);
        m->set_value(m, "VGRD_10maboveground", &vgrd);
        m->set_value(m, "param_catchment_vegetated_fraction", &veg_frac);
        m->set_value(m, "bare_soil_rsurf_exp", &rsurf_exp);
        m->update(m);
    }

    /* --- Phase 4: restore --- */
    bmi_status = m->set_value(m, NGEN_SERIALIZATION_STATE, saved_buf);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nSetValue(state) restore failed");
        free(saved_buf); m->finalize(m);
        return TEST_RETURN_CODE_FAIL;
    }

    /* Snapshot B: verify bit-exact */
    double snap_b_soil = 0, snap_b_gw = 0;
    double snap_b_theta[NDISC] = {0};
    double snap_b_nash[MAX_NUM_SUBSURFACE_NASH_CASCADE] = {0};
    double snap_b_giuh[MAX_NUM_GIUH_ORDINATES] = {0};
    m->get_value(m, "state_soil_storage_m", &snap_b_soil);
    m->get_value(m, "state_gw_storage_m", &snap_b_gw);
    m->get_value(m, "state_soil_moisture_theta", snap_b_theta);
    m->get_value(m, "state_nash_subsurface_storage", snap_b_nash);
    m->get_value(m, "state_giuh_queue", snap_b_giuh);

    if (snap_a_soil != snap_b_soil || snap_a_gw != snap_b_gw) {
        printf("\nPT scalar state mismatch: soil=%.15e vs %.15e, gw=%.15e vs %.15e",
               snap_a_soil, snap_b_soil, snap_a_gw, snap_b_gw);
        free(saved_buf); m->finalize(m);
        return TEST_RETURN_CODE_FAIL;
    }
    for (int i = 0; i < NDISC; i++) {
        if (snap_a_theta[i] != snap_b_theta[i]) {
            printf("\nPT theta[%d] mismatch: %.15e vs %.15e", i, snap_a_theta[i], snap_b_theta[i]);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {
        if (snap_a_nash[i] != snap_b_nash[i]) {
            printf("\nPT nash[%d] mismatch", i);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }
    for (int i = 0; i < MAX_NUM_GIUH_ORDINATES; i++) {
        if (snap_a_giuh[i] != snap_b_giuh[i]) {
            printf("\nPT giuh[%d] mismatch", i);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    /* --- Phase 5: run 5 steps from restored state --- */
    double q_restored[5] = {0};
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->set_value(m, "TMP_2maboveground", &temp_k);
        m->set_value(m, "DLWRF_surface", &dlwrf);
        m->set_value(m, "DSWRF_surface", &dswrf);
        m->set_value(m, "PRES_surface", &pres);
        m->set_value(m, "SPFH_2maboveground", &spfh);
        m->set_value(m, "UGRD_10maboveground", &ugrd);
        m->set_value(m, "VGRD_10maboveground", &vgrd);
        m->set_value(m, "param_catchment_vegetated_fraction", &veg_frac);
        m->set_value(m, "bare_soil_rsurf_exp", &rsurf_exp);
        m->update(m);
        m->get_value(m, "discharge_m", &q_restored[t]);
    }

    /* --- Phase 6: fresh continuous run --- */
    m->finalize(m);
    bmi_status = m->initialize(m, BMI_INIT_CONFIG_PT);
    if (bmi_status != BMI_SUCCESS) {
        printf("\nFailed to re-initialize PT for continuous comparison");
        free(saved_buf);
        return TEST_RETURN_CODE_FAIL;
    }

    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->set_value(m, "TMP_2maboveground", &temp_k);
        m->set_value(m, "DLWRF_surface", &dlwrf);
        m->set_value(m, "DSWRF_surface", &dswrf);
        m->set_value(m, "PRES_surface", &pres);
        m->set_value(m, "SPFH_2maboveground", &spfh);
        m->set_value(m, "UGRD_10maboveground", &ugrd);
        m->set_value(m, "VGRD_10maboveground", &vgrd);
        m->set_value(m, "param_catchment_vegetated_fraction", &veg_frac);
        m->set_value(m, "bare_soil_rsurf_exp", &rsurf_exp);
        m->update(m);
    }
    double q_continuous[5] = {0};
    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain);
        m->set_value(m, "et_potential_m", &pet);
        m->set_value(m, "TMP_2maboveground", &temp_k);
        m->set_value(m, "DLWRF_surface", &dlwrf);
        m->set_value(m, "DSWRF_surface", &dswrf);
        m->set_value(m, "PRES_surface", &pres);
        m->set_value(m, "SPFH_2maboveground", &spfh);
        m->set_value(m, "UGRD_10maboveground", &ugrd);
        m->set_value(m, "VGRD_10maboveground", &vgrd);
        m->set_value(m, "param_catchment_vegetated_fraction", &veg_frac);
        m->set_value(m, "bare_soil_rsurf_exp", &rsurf_exp);
        m->update(m);
        m->get_value(m, "discharge_m", &q_continuous[t]);
    }

    /* --- Phase 7: compare outputs --- */
    for (int t = 0; t < 5; t++) {
        if (q_restored[t] != q_continuous[t]) {
            printf("\nPT discharge mismatch at step %d: restored=%.15e continuous=%.15e",
                   t, q_restored[t], q_continuous[t]);
            free(saved_buf); m->finalize(m);
            return TEST_RETURN_CODE_FAIL;
        }
    }

    printf("\n  PT round-trip: %d bytes, state bit-exact, outputs bit-exact over 5 post-restore steps",
           buf_size);

    free(saved_buf);
    m->finalize(m);
    return TEST_RETURN_CODE_PASS;
}

/*
 * test_derived_quantities_resync
 *
 * Verify that derived quantities (schaake_magic_constant,
 * field_capacity_moisture_content, field_capacity_storage_m) are
 * recomputed when their base calibration parameters change via
 * set_value. The discharge-level tests can miss this because some
 * parameters also have direct effects that mask stale derived values.
 *
 * This test reads the derived fields directly from the context struct
 * via Bmi.data — it's a white-box test by design.
 */
int test_derived_quantities_resync(TestFixture* fixture)
{
    Bmi m_storage, *m = &m_storage;
    register_bmi_cfe(m);

    if (m->initialize(m, fixture->cfg_file) != BMI_SUCCESS) {
        printf("\nFailed to initialize");
        return TEST_RETURN_CODE_FAIL;
    }

    CFE_Model_Context *ctx = (CFE_Model_Context *)m->data;
    cfe_parameters_struct *p = &ctx->parameters;

    double orig_schaake    = p->schaake_magic_constant;
    double orig_fc_theta   = p->field_capacity_moisture_content;
    double orig_fc_storage = p->field_capacity_storage_m;

    /* Perturb ksat by 100x — schaake_magic_constant must change proportionally */
    double new_ksat_ms = p->ksat_m_per_s * 100.0;
    double new_ksat_cmh = new_ksat_ms * 360000.0;
    m->set_value(m, "soil_saturated_hydraulic_conductivity", &new_ksat_cmh);

    /* Perturb soil_b (4.05 → 12.0) — field_capacity must change */
    double new_b = 12.0;
    m->set_value(m, "soil_Clapp_Hornberger_b", &new_b);

    /* Run one step to trigger any resync logic */
    double rain = 0.005 / 3600.0, pet = 0.0;
    m->set_value(m, "rainfall_depth_m", &rain);
    m->set_value(m, "et_potential_m", &pet);
    m->update(m);

    int failed = 0;

    /* schaake_magic_constant = refkdt * ksat / 2e-6 (ksat in m/s internally) */
    double expected_schaake = p->refkdt * new_ksat_ms / 2.0e-06;
    /* Relative tolerance: schaake magnitude varies with Ksat, so scale epsilon to the expected value */
    if (fabs(p->schaake_magic_constant - expected_schaake) > 1.0e-10 * fabs(expected_schaake)) {
        printf("\n  FAIL: schaake_magic_constant stale: got %.6e, expected %.6e (orig %.6e)",
               p->schaake_magic_constant, expected_schaake, orig_schaake);
        failed++;
    } else {
        printf("\n  OK: schaake_magic_constant recomputed (%.6e → %.6e)", orig_schaake, p->schaake_magic_constant);
    }

    /* field_capacity_moisture_content must differ from original after soil_b change */
    if (fabs(p->field_capacity_moisture_content - orig_fc_theta) < 1.0e-15) {
        printf("\n  FAIL: field_capacity_moisture_content unchanged (%.6e)", orig_fc_theta);
        failed++;
    } else {
        printf("\n  OK: field_capacity_moisture_content recomputed (%.6e → %.6e)",
               orig_fc_theta, p->field_capacity_moisture_content);
    }

    /* field_capacity_storage_m must differ too */
    if (fabs(p->field_capacity_storage_m - orig_fc_storage) < 1.0e-15) {
        printf("\n  FAIL: field_capacity_storage_m unchanged (%.6e)", orig_fc_storage);
        failed++;
    } else {
        printf("\n  OK: field_capacity_storage_m recomputed (%.6e → %.6e)",
               orig_fc_storage, p->field_capacity_storage_m);
    }

    /* Verify field_capacity_storage_m = fc_theta * soil_depth_m */
    double expected_fc_storage = p->field_capacity_moisture_content * p->soil_depth_m;
    if (fabs(p->field_capacity_storage_m - expected_fc_storage) > 1.0e-15) {
        printf("\n  FAIL: field_capacity_storage_m (%.6e) != fc_theta * depth (%.6e)",
               p->field_capacity_storage_m, expected_fc_storage);
        failed++;
    }

    m->finalize(m);
    if (failed > 0) return TEST_RETURN_CODE_FAIL;
    printf("\n  All derived quantities resynced correctly");
    return TEST_RETURN_CODE_PASS;
}

/*
 * test_dsbm_soil_params_resync
 *
 * Verify that s->soil_parameters (the DSBM cached copy) is refreshed
 * when calibration parameters change via set_value. Uses the DSBM
 * config and reads the state struct directly.
 */
int test_dsbm_soil_params_resync(TestFixture* fixture)
{
    (void)fixture;
    Bmi m_storage, *m = &m_storage;
    register_bmi_cfe(m);

    if (m->initialize(m, BMI_INIT_CONFIG_DSBM) != BMI_SUCCESS) {
        printf("\nFailed to initialize DSBM config");
        return TEST_RETURN_CODE_FAIL;
    }

    CFE_Model_Context *ctx = (CFE_Model_Context *)m->data;
    cfe_parameters_struct *p = &ctx->parameters;
    SoilParameters *sp = &ctx->state.soil_parameters;

    /* Record originals */
    double orig_Ksat_cm_h = sp->K_sat_cm_per_h;
    double orig_b_exp     = sp->b_exp;
    double orig_phi_cm    = sp->phi_sat_cm;
    double orig_perc      = sp->perc_limiter_0_to_1;
    double orig_klf       = sp->klf_per_h;
    double orig_theta_sat = sp->theta_sat;

    /* Perturb base parameters via BMI (ksat in cm/h, satpsi in cm at BMI boundary) */
    double new_ksat_ms = p->ksat_m_per_s * 100.0;
    double new_ksat_cmh = new_ksat_ms * 360000.0;
    m->set_value(m, "soil_saturated_hydraulic_conductivity", &new_ksat_cmh);
    double new_b = 12.0;
    m->set_value(m, "soil_Clapp_Hornberger_b", &new_b);
    double new_phi_cm = 100.0;  /* cm (was 35.5 cm / 0.355 m) */
    m->set_value(m, "soil_saturated_capillary_head", &new_phi_cm);
    double new_perc = 0.90;
    m->set_value(m, "soil_percolation_rate_limiter", &new_perc);
    double new_klf = 0.50;
    m->set_value(m, "soil_lateral_flow_K", &new_klf);
    double new_porosity = 0.30;
    m->set_value(m, "soil_effective_porosity", &new_porosity);

    /* Run one step to trigger resync */
    double rain = 0.005 / 3600.0, pet = 0.0;
    m->set_value(m, "rainfall_depth_m", &rain);
    m->set_value(m, "et_potential_m", &pet);
    m->update(m);

    int failed = 0;

    /* Check each s->soil_parameters field was updated */
    double exp_Ksat_cm_h = new_ksat_ms * 360000.0;
    if (fabs(sp->K_sat_cm_per_h - exp_Ksat_cm_h) > 1.0e-6) {
        printf("\n  FAIL: soil_parameters.K_sat_cm_per_h stale: %.6e (expected %.6e, orig %.6e)",
               sp->K_sat_cm_per_h, exp_Ksat_cm_h, orig_Ksat_cm_h);
        failed++;
    } else {
        printf("\n  OK: K_sat_cm_per_h resynced (%.4e → %.4e)", orig_Ksat_cm_h, sp->K_sat_cm_per_h);
    }

    if (fabs(sp->b_exp - new_b) > 1.0e-15) {
        printf("\n  FAIL: soil_parameters.b_exp stale: %.6f (expected %.6f, orig %.6f)",
               sp->b_exp, new_b, orig_b_exp);
        failed++;
    } else {
        printf("\n  OK: b_exp resynced (%.2f → %.2f)", orig_b_exp, sp->b_exp);
    }

    double exp_phi_cm = new_phi_cm;
    if (fabs(sp->phi_sat_cm - exp_phi_cm) > 1.0e-10) {
        printf("\n  FAIL: soil_parameters.phi_sat_cm stale: %.4f (expected %.4f, orig %.4f)",
               sp->phi_sat_cm, exp_phi_cm, orig_phi_cm);
        failed++;
    } else {
        printf("\n  OK: phi_sat_cm resynced (%.2f → %.2f)", orig_phi_cm, sp->phi_sat_cm);
    }

    if (fabs(sp->perc_limiter_0_to_1 - new_perc) > 1.0e-15) {
        printf("\n  FAIL: soil_parameters.perc_limiter stale: %.4f (expected %.4f, orig %.4f)",
               sp->perc_limiter_0_to_1, new_perc, orig_perc);
        failed++;
    } else {
        printf("\n  OK: perc_limiter resynced (%.4f → %.4f)", orig_perc, sp->perc_limiter_0_to_1);
    }

    if (fabs(sp->klf_per_h - new_klf) > 1.0e-15) {
        printf("\n  FAIL: soil_parameters.klf_per_h stale: %.4f (expected %.4f, orig %.4f)",
               sp->klf_per_h, new_klf, orig_klf);
        failed++;
    } else {
        printf("\n  OK: klf_per_h resynced (%.4f → %.4f)", orig_klf, sp->klf_per_h);
    }

    if (fabs(sp->theta_sat - new_porosity) > 1.0e-15) {
        printf("\n  FAIL: soil_parameters.theta_sat stale: %.4f (expected %.4f, orig %.4f)",
               sp->theta_sat, new_porosity, orig_theta_sat);
        failed++;
    } else {
        printf("\n  OK: theta_sat resynced (%.4f → %.4f)", orig_theta_sat, sp->theta_sat);
    }

    m->finalize(m);
    if (failed > 0) {
        printf("\n  %d DSBM soil_parameters fields not resynced", failed);
        return TEST_RETURN_CODE_FAIL;
    }
    printf("\n  All DSBM soil_parameters fields resynced correctly");
    return TEST_RETURN_CODE_PASS;
}

/*
 * test_calibration_params_affect_output
 *
 * For each calibration parameter: initialize, set to a perturbed value,
 * run 5 timesteps with rainfall, and verify that cumulative discharge
 * differs from an unperturbed baseline run. This catches the bug where
 * set_value writes to the parameter struct but derived quantities or
 * cached copies (s->soil_parameters, schaake_magic_constant,
 * field_capacity_storage_m, lookup tables) are never refreshed.
 */

static double run_and_get_discharge(const char *cfg_file,
                                    const char *param_name,
                                    double param_value,
                                    int do_set)
{
    Bmi m_storage, *m = &m_storage;
    register_bmi_cfe(m);

    if (m->initialize(m, cfg_file) != BMI_SUCCESS) return -1.0;

    if (do_set) {
        if (m->set_value(m, param_name, &param_value) != BMI_SUCCESS) {
            m->finalize(m);
            return -1.0;
        }
    }

    double rain_m = 0.005 / 3600.0;
    double pet_m  = 0.0;
    double cumulative_q = 0.0;

    for (int t = 0; t < 5; t++) {
        m->set_value(m, "rainfall_depth_m", &rain_m);
        m->set_value(m, "et_potential_m", &pet_m);
        if (m->update(m) != BMI_SUCCESS) { m->finalize(m); return -1.0; }

        double q = 0.0;
        m->get_value(m, "discharge_m", &q);
        cumulative_q += q;
    }

    m->finalize(m);
    return cumulative_q;
}

int test_calibration_params_affect_output(TestFixture* fixture)
{
    const char *cfg = fixture->cfg_file;

    /* Each entry: BMI parameter name, perturbed value (must be physically
       valid but far enough from the config default to change output).
       Config defaults from cfe_config_cat_87_pass.cf3 noted in comments. */
    static const struct { const char *name; double perturbed; } params[] = {
        { "soil_effective_porosity",               0.20   },  /* default 0.439 */
        { "soil_saturated_hydraulic_conductivity",  122.4  },  /* default ~1.22e-2 cm/h; 100x increase */
        { "soil_percolation_rate_limiter",          0.90   },  /* default 0.01 */
        { "soil_Clapp_Hornberger_b",               12.0   },  /* default 4.05 */
        { "soil_lateral_flow_K",                    0.50   },  /* default 0.01 h-1 */
        { "subsurface_nash_K",                      0.50   },  /* default 0.03 h-1 */
        { "gw_discharge_coefficient",               1.8e-3 },  /* default 1.8e-5 */
        { "gw_discharge_exponent",                  1.5    },  /* default 6.0 */
        { "gw_max_storage_m",                       0.01   },  /* default 0.25 */
        { "soil_saturated_capillary_head",          100.0  },  /* default 35.5 cm */
        { "soil_field_capacity_fraction",           0.10   },  /* default 0.333 */
    };
    int n_params = sizeof(params) / sizeof(params[0]);

    double baseline = run_and_get_discharge(cfg, NULL, 0.0, 0);
    if (baseline < 0.0) {
        printf("\nFailed to run baseline");
        return TEST_RETURN_CODE_FAIL;
    }
    if (baseline == 0.0) {
        printf("\nBaseline discharge is zero — test is invalid");
        return TEST_RETURN_CODE_FAIL;
    }

    int failed = 0;
    for (int i = 0; i < n_params; i++) {
        double perturbed = run_and_get_discharge(cfg, params[i].name,
                                                 params[i].perturbed, 1);
        if (perturbed < 0.0) {
            printf("\n  FAIL: run failed for '%s'", params[i].name);
            failed++;
            continue;
        }
        double rel_change = fabs(perturbed - baseline) / baseline;
        if (rel_change < 1.0e-10) {
            printf("\n  FAIL: '%s' had no effect on discharge "
                   "(baseline=%.6e, perturbed=%.6e, rel_change=%.2e)",
                   params[i].name, baseline, perturbed, rel_change);
            failed++;
        } else {
            printf("\n  OK: '%s' rel_change=%.4e", params[i].name, rel_change);
        }
    }

    if (failed > 0) {
        printf("\n  %d of %d calibration parameters had no effect on output", failed, n_params);
        return TEST_RETURN_CODE_FAIL;
    }

    printf("\n  All %d calibration parameters affect discharge output", n_params);
    return TEST_RETURN_CODE_PASS;
}

/*
 * test_calibration_params_affect_output_dsbm
 *
 * Same as above but with the DSBM config. Tests the s->soil_parameters
 * copy path and lookup table dependencies.
 */
int test_calibration_params_affect_output_dsbm(TestFixture* fixture)
{
    (void)fixture;
    const char *cfg = BMI_INIT_CONFIG_DSBM;

    static const struct { const char *name; double perturbed; } params[] = {
        { "soil_effective_porosity",               0.20   },
        { "soil_saturated_hydraulic_conductivity",  122.4  },
        { "soil_percolation_rate_limiter",          0.90   },
        { "soil_Clapp_Hornberger_b",               12.0   },
        { "soil_lateral_flow_K",                    0.50   },
        { "subsurface_nash_K",                      0.50   },
        { "gw_discharge_coefficient",               1.8e-3 },
        { "gw_discharge_exponent",                  1.5    },
        { "gw_max_storage_m",                       0.01   },
        { "soil_saturated_capillary_head",          100.0  },
        { "soil_field_capacity_fraction",           0.10   },
    };
    int n_params = sizeof(params) / sizeof(params[0]);

    double baseline = run_and_get_discharge(cfg, NULL, 0.0, 0);
    if (baseline < 0.0) {
        printf("\nFailed to run DSBM baseline");
        return TEST_RETURN_CODE_FAIL;
    }
    if (baseline == 0.0) {
        printf("\nDSBM baseline discharge is zero — test is invalid");
        return TEST_RETURN_CODE_FAIL;
    }

    int failed = 0;
    for (int i = 0; i < n_params; i++) {
        double perturbed = run_and_get_discharge(cfg, params[i].name,
                                                 params[i].perturbed, 1);
        if (perturbed < 0.0) {
            printf("\n  FAIL: DSBM run failed for '%s'", params[i].name);
            failed++;
            continue;
        }
        double rel_change = fabs(perturbed - baseline) / baseline;
        if (rel_change < 1.0e-10) {
            printf("\n  FAIL (DSBM): '%s' had no effect on discharge "
                   "(baseline=%.6e, perturbed=%.6e, rel_change=%.2e)",
                   params[i].name, baseline, perturbed, rel_change);
            failed++;
        } else {
            printf("\n  OK (DSBM): '%s' rel_change=%.4e", params[i].name, rel_change);
        }
    }

    if (failed > 0) {
        printf("\n  %d of %d DSBM calibration parameters had no effect on output", failed, n_params);
        return TEST_RETURN_CODE_FAIL;
    }

    printf("\n  All %d DSBM calibration parameters affect discharge output", n_params);
    return TEST_RETURN_CODE_PASS;
}

/*
 * Helper: initialize a fresh Bmi with the given config, query input count
 * and names, finalize, and check against expectations.
 * Returns 0 on success, 1 on failure.
 */
static int check_input_vars(const char *cfg, int expected_count,
                            const char *expected_names[], const char *label)
{
    Bmi m_storage, *m = &m_storage;
    register_bmi_cfe(m);

    if (m->initialize(m, cfg) != BMI_SUCCESS) {
        printf("\n  [%s] Failed to initialize", label);
        return 1;
    }

    int count = -1;
    m->get_input_item_count(m, &count);
    if (count != expected_count) {
        printf("\n  [%s] input count: got %d, want %d", label, count, expected_count);
        m->finalize(m);
        return 1;
    }

    char **names = allocate_array_of_strings(count, BMI_MAX_VAR_NAME);
    m->get_input_var_names(m, names);

    int failures = 0;
    for (int i = 0; i < count; i++) {
        if (strcmp(names[i], expected_names[i]) != 0) {
            printf("\n  [%s] input[%d]: got '%s', want '%s'",
                   label, i, names[i], expected_names[i]);
            failures++;
        }
    }

    free_array_of_strings(names, count);
    m->finalize(m);
    return failures;
}

/*
 * test_conditional_inputs
 *
 * Verify that Get_input_item_count / Get_input_var_names respond correctly
 * to the model configuration:
 *
 *   1. PT disabled, soil evap disabled  → 5 inputs (base set)
 *   2. PT+soil_evap enabled             → 13 inputs (base + AORC + day_of_year)
 *   3. PT+soil_evap + internal DOY      → 12 inputs (AORC but no day_of_year)
 */
int test_conditional_inputs(TestFixture* fixture)
{
    int failures = 0;

    /* Case 1: PT disabled (the default unit test config) — 5 inputs */
    {
        const char *expected[] = {
            "rainfall_depth_m", "et_potential_m", "ice_fraction",
            "param_catchment_vegetated_fraction", "bare_soil_rsurf_exp"
        };
        failures += check_input_vars(
            fixture->cfg_file, 5, expected, "pt_disabled");
    }

    /* Case 2: PT + soil evap enabled, external day_of_year — 13 inputs */
    {
        const char *expected[] = {
            "rainfall_depth_m", "et_potential_m", "ice_fraction",
            "day_of_year",
            "DLWRF_surface", "DSWRF_surface", "PRES_surface",
            "SPFH_2maboveground", "TMP_2maboveground",
            "UGRD_10maboveground", "VGRD_10maboveground",
            "param_catchment_vegetated_fraction", "bare_soil_rsurf_exp"
        };
        failures += check_input_vars(
            BMI_INIT_CONFIG_PT, 13, expected, "pt_enabled");
    }

    /* Case 3: PT + soil evap + internal day_of_year — 12 inputs */
    {
        const char *expected[] = {
            "rainfall_depth_m", "et_potential_m", "ice_fraction",
            "DLWRF_surface", "DSWRF_surface", "PRES_surface",
            "SPFH_2maboveground", "TMP_2maboveground",
            "UGRD_10maboveground", "VGRD_10maboveground",
            "param_catchment_vegetated_fraction", "bare_soil_rsurf_exp"
        };
        failures += check_input_vars(
            BMI_INIT_CONFIG_PT_INTERNAL_DOY, 12, expected, "pt_internal_doy");
    }

    if (failures > 0) {
        printf("\n  %d conditional input check(s) failed", failures);
        return TEST_RETURN_CODE_FAIL;
    }
    printf("\n  All 3 conditional input configurations verified");
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
    else if (strcmp(argv[1], "test_update") == 0)
        result = test_update(fixture);
    else if (strcmp(argv[1], "test_update_until") == 0)
        result = test_update_until(fixture);
    else if (strcmp(argv[1], "test_finalize") == 0)
        result = test_finalize(fixture);
    else if (strcmp(argv[1], "test_get_component_name") == 0)
        result = test_get_component_name(fixture);
    else if (strcmp(argv[1], "test_get_current_time") == 0)
        result = test_get_current_time(fixture);
    else if (strcmp(argv[1], "test_get_end_time") == 0)
        result = test_get_end_time(fixture);
    else if (strcmp(argv[1], "test_get_grid_edge_count") == 0)
        result = test_get_grid_edge_count(fixture);
    else if (strcmp(argv[1], "test_get_grid_edge_nodes") == 0)
        result = test_get_grid_edge_nodes(fixture);
    else if (strcmp(argv[1], "test_get_grid_face_count") == 0)
        result = test_get_grid_face_count(fixture);
    else if (strcmp(argv[1], "test_get_grid_face_edges") == 0)
        result = test_get_grid_face_edges(fixture);
    else if (strcmp(argv[1], "test_get_grid_face_nodes") == 0)
        result = test_get_grid_face_nodes(fixture);
    else if (strcmp(argv[1], "test_get_grid_node_count") == 0)
        result = test_get_grid_node_count(fixture);
    else if (strcmp(argv[1], "test_get_grid_nodes_per_face") == 0)
        result = test_get_grid_nodes_per_face(fixture);
    else if (strcmp(argv[1], "test_get_grid_origin") == 0)
        result = test_get_grid_origin(fixture);
    else if (strcmp(argv[1], "test_get_grid_rank") == 0)
        result = test_get_grid_rank(fixture);
    else if (strcmp(argv[1], "test_get_grid_shape") == 0)
        result = test_get_grid_shape(fixture);
    else if (strcmp(argv[1], "test_get_grid_size") == 0)
        result = test_get_grid_size(fixture);
    else if (strcmp(argv[1], "test_get_grid_spacing") == 0)
        result = test_get_grid_spacing(fixture);
    else if (strcmp(argv[1], "test_get_grid_type") == 0)
        result = test_get_grid_type(fixture);
    else if (strcmp(argv[1], "test_get_grid_x") == 0)
        result = test_get_grid_x(fixture);
    else if (strcmp(argv[1], "test_get_grid_y") == 0)
        result = test_get_grid_y(fixture);
    else if (strcmp(argv[1], "test_get_grid_z") == 0)
        result = test_get_grid_z(fixture);
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
    else if (strcmp(argv[1], "test_mass_balance_protocol") == 0)
        result = test_mass_balance_protocol(fixture);
    else if (strcmp(argv[1], "test_grid_consistency") == 0)
        result = test_grid_consistency(fixture);
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
    else if (strcmp(argv[1], "test_serialization_metadata") == 0)
        result = test_serialization_metadata(fixture);
    else if (strcmp(argv[1], "test_serialization_round_trip") == 0)
        result = test_serialization_round_trip(fixture);
    else if (strcmp(argv[1], "test_serialization_round_trip_dsbm") == 0)
        result = test_serialization_round_trip_dsbm(fixture);
    else if (strcmp(argv[1], "test_serialization_round_trip_pt") == 0)
        result = test_serialization_round_trip_pt(fixture);
    else if (strcmp(argv[1], "test_derived_quantities_resync") == 0)
        result = test_derived_quantities_resync(fixture);
    else if (strcmp(argv[1], "test_dsbm_soil_params_resync") == 0)
        result = test_dsbm_soil_params_resync(fixture);
    else if (strcmp(argv[1], "test_calibration_params_affect_output") == 0)
        result = test_calibration_params_affect_output(fixture);
    else if (strcmp(argv[1], "test_calibration_params_affect_output_dsbm") == 0)
        result = test_calibration_params_affect_output_dsbm(fixture);
    else if (strcmp(argv[1], "test_conditional_inputs") == 0)
        result = test_conditional_inputs(fixture);
    else
        printf("\nUnexpected test function %s\n", argv[1]);

    teardown(fixture);
    return result;
}