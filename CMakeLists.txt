cmake_minimum_required(VERSION 3.10)

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red         "${Esc}[32m")

# Build options
# NGEN=ON (default): build shared library for ngen framework integration
# STANDALONE=ON: additionally build the standalone non-BMI driver executable
option(NGEN "Build CFE BMI shared library (ON by default)" ON)
option(STANDALONE "Build standalone (non-BMI) CFE driver executable" OFF)

# A numeric verbosity level can be set for debug output printed/logged during execution
# Use 1 by default when build type is "Debug"
if (NOT DEFINED DEBUG_VERBOSITY AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(DEBUG_VERBOSITY 1)
# Use 0 by default for any build types other than "Debug"
elseif (NOT DEFINED DEBUG_VERBOSITY)
    set(DEBUG_VERBOSITY 0)
# Catch any "ON" or equivalent, excluding all integers other than "1" (which we want to handle separately)
elseif (${DEBUG_VERBOSITY} AND NOT ${DEBUG_VERBOSITY} GREATER "1" AND NOT ${DEBUG_VERBOSITY} LESS_EQUAL "0")
    set(DEBUG_VERBOSITY 1)
# Should catch any other values except an integer greater than 1 (which we want to use directly)
elseif (NOT ${DEBUG_VERBOSITY} GREATER "1")
    set(DEBUG_VERBOSITY 0)
endif ()
# While not an explicit condition, DEBUG_VERBOSITY is left alone if set to some integer greater than 1

# TODO: this should ideally be applied to a specific target, but how targets are defined and set up probably needs to be adjusted first
add_compile_definitions(CFE_DEBUG=${DEBUG_VERBOSITY})

# set the project name
project(cfebmi VERSION 3.0.0 DESCRIPTION "OWP CFE BMI Module Shared Library")

IF (CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF (CMAKE_BUILD_TYPE MATCHES Debug)

message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# CFE shared library

set(CFE_CORE_SOURCES
    ./src/bmi_cfe.c
    ./src/cfe_serialize.c
    ./src/cfe.c
    ./src/giuh.c
    ./src/conceptual_reservoir.c
    ./src/nash_cascade.c
    ./src/parser_helpers.c
    ./src/cfe_helpers.c
    ./src/cfe_context.c
    ./src/cfe_driver_utils.c
    ./src/soil_helpers.c
    ./src/discrete_soil_moisture_stateless.c
    ./src/cfe_pet_priestley_taylor.c
    ./src/cfe_soil_skin_temperature.c
    ./src/calculate_bare_soil_evap.c
)

if (WIN32)
    add_library(cfebmi ${CFE_CORE_SOURCES})
else ()
    add_library(cfebmi SHARED ${CFE_CORE_SOURCES})
endif ()

# On some platforms, the linker needs to be explicitly told to search for defs for standard C math functions
target_link_libraries(cfebmi PRIVATE m)

target_include_directories(cfebmi PRIVATE include)

set_target_properties(cfebmi PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(cfebmi PROPERTIES PUBLIC_HEADER ./include/bmi_cfe.h)

####################################################################################################################
# Config migration utility (v2 -> v3, self-contained, no library dependency)
add_executable(cfe_migrate_config ./tools/cfe_migrate_config.c)
target_link_libraries(cfe_migrate_config PRIVATE m)

####################################################################################################################
# CFE BMI driver executable (used by integration tests and standalone operation)
add_executable(cfe_bmi_driver ./src/cfe_bmi_driver.c)
target_include_directories(cfe_bmi_driver PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(cfe_bmi_driver cfebmi)

####################################################################################################################
# Standalone (non-BMI) driver executable — always built (used by restart integration tests)
add_executable(cfe_main_driver ./src/cfe_main_driver.c)
target_include_directories(cfe_main_driver PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(cfe_main_driver cfebmi)

####################################################################################################################
# Testing
enable_testing()

####################################################################################################################
# Include executable and (single) test for the legacy batch test file
add_executable(combined_bmi_funcs_test test/combined_bmi_funcs_test.c)
target_include_directories(combined_bmi_funcs_test PUBLIC ${PROJECT_SOURCE_DIR}/include test/include)
target_link_libraries(combined_bmi_funcs_test cfebmi)
add_test(NAME combined_bmi_funcs_test_1 COMMAND combined_bmi_funcs_test ${PROJECT_SOURCE_DIR}/configs/cfe_config_cat_87_pass.cf3)
# End combined_bmi_funcs_test tests
####################################################################################################################

####################################################################################################################
# Include executable and individual, isolated unit tests for BMI related functions
# Control which test is run via custom implemented arg parsing within the executable
add_executable(test_bmi_model ./test/test_bmi_model.c test/bmi_test_utils.c test/bmi_test_utils.c test/general_test_utils.c)
target_include_directories(test_bmi_model PUBLIC ${PROJECT_SOURCE_DIR}/include test/include)
target_link_libraries(test_bmi_model cfebmi)
target_compile_definitions(test_bmi_model PUBLIC
    BMI_INIT_CONFIG_EX_1="${PROJECT_SOURCE_DIR}/configs/cfe_config_cat_87_pass.cf3"
    BMI_INIT_CONFIG_DSBM="${PROJECT_SOURCE_DIR}/configs/clean_config.cf3"
    BMI_INIT_CONFIG_PT="${PROJECT_SOURCE_DIR}/configs/test_pt_enabled.cf3"
    BMI_INIT_CONFIG_PT_INTERNAL_DOY="${PROJECT_SOURCE_DIR}/configs/test_pt_internal_doy.cf3")

# For the 'test_bmi_model' executable, the test example can be explicitly passed (although for now there is only 1)
add_test(NAME test_initialize_1 COMMAND test_bmi_model test_initialize 1)
# Note that, when using the first test example, the value be left off/implied for convenience
add_test(NAME test_initialize_1_a COMMAND test_bmi_model test_initialize)
add_test(NAME test_update_1 COMMAND test_bmi_model test_update)
add_test(NAME test_update_until_1 COMMAND test_bmi_model test_update_until)
add_test(NAME test_finalize_1 COMMAND test_bmi_model test_finalize)

add_test(NAME test_get_component_name_1 COMMAND test_bmi_model test_get_component_name)
add_test(NAME test_get_input_item_count_1 COMMAND test_bmi_model test_get_input_item_count)
add_test(NAME test_get_input_var_names_1 COMMAND test_bmi_model test_get_input_var_names)
add_test(NAME test_get_output_item_count_1 COMMAND test_bmi_model test_get_output_item_count)
add_test(NAME test_get_output_var_names_1 COMMAND test_bmi_model test_get_output_var_names)
add_test(NAME test_get_var_grid_1 COMMAND test_bmi_model test_get_var_grid)
add_test(NAME test_get_var_itemsize_1 COMMAND test_bmi_model test_get_var_itemsize)
add_test(NAME test_get_var_location_1 COMMAND test_bmi_model test_get_var_location)
add_test(NAME test_get_var_units_1 COMMAND test_bmi_model test_get_var_units)
add_test(NAME test_get_var_type_1 COMMAND test_bmi_model test_get_var_type)
# TODO: While others could also benefit, get_var_nbytes in particular needs a second test example that covers the
# TODO:     other case for the "soil_moisture_profile" variable (see comments in test_get_var_nbytes function)
add_test(NAME test_get_var_nbytes_1 COMMAND test_bmi_model test_get_var_nbytes)
add_test(NAME test_get_grid_rank_1 COMMAND test_bmi_model test_get_grid_rank)
add_test(NAME test_get_grid_size_1 COMMAND test_bmi_model test_get_grid_size)
add_test(NAME test_get_grid_type_1 COMMAND test_bmi_model test_get_grid_type)
add_test(NAME test_get_start_time_1 COMMAND test_bmi_model test_get_start_time)
# TODO: (later) get_end_time probably needs another test scenario when the module is handling it's own forcings
add_test(NAME test_get_end_time_1 COMMAND test_bmi_model test_get_end_time)
add_test(NAME test_get_time_step_1 COMMAND test_bmi_model test_get_time_step)
add_test(NAME test_get_time_units_1 COMMAND test_bmi_model test_get_time_units)
add_test(NAME test_get_current_time_1 COMMAND test_bmi_model test_get_current_time)
add_test(NAME test_set_value_1 COMMAND test_bmi_model test_set_value)
add_test(NAME test_set_value_at_indices_1 COMMAND test_bmi_model test_set_value_at_indices)
add_test(NAME test_get_value_1 COMMAND test_bmi_model test_get_value)
add_test(NAME test_get_value_at_indices_1 COMMAND test_bmi_model test_get_value_at_indices)
add_test(NAME test_get_value_ptr_1 COMMAND test_bmi_model test_get_value_ptr)
add_test(NAME test_mass_balance_protocol_1 COMMAND test_bmi_model test_mass_balance_protocol)
add_test(NAME test_grid_consistency_1 COMMAND test_bmi_model test_grid_consistency)

# These are effectively not implemented, but they should strictly speaking be implemented and return a failure code
add_test(NAME test_get_grid_edge_count_1 COMMAND test_bmi_model test_get_grid_edge_count)
add_test(NAME test_get_grid_edge_nodes_1 COMMAND test_bmi_model test_get_grid_edge_nodes)
add_test(NAME test_get_grid_face_count_1 COMMAND test_bmi_model test_get_grid_face_count)
add_test(NAME test_get_grid_face_edges_1 COMMAND test_bmi_model test_get_grid_face_edges)
add_test(NAME test_get_grid_face_nodes_1 COMMAND test_bmi_model test_get_grid_face_nodes)
add_test(NAME test_get_grid_node_count_1 COMMAND test_bmi_model test_get_grid_node_count)
add_test(NAME test_get_grid_nodes_per_face_1 COMMAND test_bmi_model test_get_grid_nodes_per_face)
add_test(NAME test_get_grid_origin_1 COMMAND test_bmi_model test_get_grid_origin)
add_test(NAME test_get_grid_shape_1 COMMAND test_bmi_model test_get_grid_shape)
add_test(NAME test_get_grid_spacing_1 COMMAND test_bmi_model test_get_grid_spacing)
add_test(NAME test_get_grid_x_1 COMMAND test_bmi_model test_get_grid_x)
add_test(NAME test_get_grid_y_1 COMMAND test_bmi_model test_get_grid_y)
add_test(NAME test_get_grid_z_1 COMMAND test_bmi_model test_get_grid_z)
add_test(NAME test_serialization_metadata_1 COMMAND test_bmi_model test_serialization_metadata)
add_test(NAME test_serialization_round_trip_1 COMMAND test_bmi_model test_serialization_round_trip)
add_test(NAME test_serialization_round_trip_dsbm_1 COMMAND test_bmi_model test_serialization_round_trip_dsbm)
add_test(NAME test_serialization_round_trip_pt_1 COMMAND test_bmi_model test_serialization_round_trip_pt)
add_test(NAME test_derived_quantities_resync COMMAND test_bmi_model test_derived_quantities_resync)
add_test(NAME test_dsbm_soil_params_resync COMMAND test_bmi_model test_dsbm_soil_params_resync)
add_test(NAME test_calibration_params_affect_output COMMAND test_bmi_model test_calibration_params_affect_output)
add_test(NAME test_calibration_params_affect_output_dsbm COMMAND test_bmi_model test_calibration_params_affect_output_dsbm)
add_test(NAME test_conditional_inputs COMMAND test_bmi_model test_conditional_inputs)
# End test_bmi_model tests
####################################################################################################################

####################################################################################################################
# Serialization unit tests (test the serialize/deserialize layer directly, not through BMI)
add_executable(test_serialization ./test/test_serialization.c)
target_include_directories(test_serialization PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(test_serialization cfebmi)
target_compile_definitions(test_serialization PUBLIC
    TEST_SERIALIZATION_CONFIG="${PROJECT_SOURCE_DIR}/configs/cfe_config_cat_87_pass.cf3")

add_test(NAME test_serialize_create_1 COMMAND test_serialization test_serialize_create)
add_test(NAME test_deserialize_failures_1 COMMAND test_serialization test_deserialize_failures)
# End test_serialization tests
####################################################################################################################

####################################################################################################################
# Integration tests: run driver with golden configs and compare numerical output
set(INTEGRATION_TEST_SCRIPT ${PROJECT_SOURCE_DIR}/test/run_integration_test.sh)
set(FORCING_FILE ${PROJECT_SOURCE_DIR}/forcings/cat87_01Dec2015.csv)

add_test(NAME integration_v2_migrated
    COMMAND bash ${INTEGRATION_TEST_SCRIPT}
        $<TARGET_FILE:cfe_bmi_driver>
        ${PROJECT_SOURCE_DIR}/configs/migrated_legacy_cat87.cf3
        ${FORCING_FILE}
        ${PROJECT_SOURCE_DIR}/test/golden/v2
        1e-10
)
set_tests_properties(integration_v2_migrated PROPERTIES LABELS "integration")

add_test(NAME integration_v3_dsbm
    COMMAND bash ${INTEGRATION_TEST_SCRIPT}
        $<TARGET_FILE:cfe_bmi_driver>
        ${PROJECT_SOURCE_DIR}/configs/bmi_config_cat87_v3.cf3
        ${FORCING_FILE}
        ${PROJECT_SOURCE_DIR}/test/golden/v3
        1e-10
)
set_tests_properties(integration_v3_dsbm PROPERTIES LABELS "integration")

add_test(NAME integration_v3_dry_start
    COMMAND bash ${INTEGRATION_TEST_SCRIPT}
        $<TARGET_FILE:cfe_bmi_driver>
        ${PROJECT_SOURCE_DIR}/configs/laramie_dry_start_v3.cf3
        ${PROJECT_SOURCE_DIR}/forcings/Laramie_dry_start_30day.csv
        ${PROJECT_SOURCE_DIR}/test/golden/v3_dry_start
        1e-10
)
set_tests_properties(integration_v3_dry_start PROPERTIES LABELS "integration")

# Restart test: split forcing, run two halves via hotstart, compare to continuous run
set(RESTART_TEST_SCRIPT ${PROJECT_SOURCE_DIR}/test/run_restart_test.sh)

add_test(NAME restart_standalone_v3
    COMMAND bash ${RESTART_TEST_SCRIPT}
        $<TARGET_FILE:cfe_main_driver>
        ${PROJECT_SOURCE_DIR}/configs/clean_config.cf3
        ${FORCING_FILE}
        360
        1e-10
)
set_tests_properties(restart_standalone_v3 PROPERTIES LABELS "integration")

add_test(NAME restart_standalone_v3_pt
    COMMAND bash ${RESTART_TEST_SCRIPT}
        $<TARGET_FILE:cfe_main_driver>
        ${PROJECT_SOURCE_DIR}/configs/test_pt_enabled.cf3
        ${FORCING_FILE}
        360
        1e-10
)
set_tests_properties(restart_standalone_v3_pt PROPERTIES LABELS "integration")
# End integration tests
####################################################################################################################

include(GNUInstallDirs)

install(TARGETS cfebmi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(cfebmi.pc.in cfebmi.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/cfebmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
