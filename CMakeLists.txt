cmake_minimum_required(VERSION 3.10)

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red         "${Esc}[32m")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# module setup options

option(BASE "BASE" OFF)
option(FORCING "FORCING" OFF)
option(FORCINGPET "FORCINGPET" OFF)
option(AETROOTZONE "AETROOTZONE" OFF)
option(NGEN "NGEN" OFF)
option(UNITTEST "UNITTEST" OFF)

if ((NOT BASE) AND (NOT FORCING) AND (NOT FORCINGPET) AND (NOT AETROOTZONE) AND (NOT NGEN) AND (NOT UNITTEST))
    message("${Red}Options: BASE, FORCING, FORCINGPET, AETROOTZONE, NGEN, UNITTEST" ${ColourReset})
    message(FATAL_ERROR "Invalid option is provided, CMake will exit.")
endif ()

if (BASE)
    message("${Red} CFE BUILD 'BASE' CASE!${ColourReset}")
    set(exe_name "cfe_base")
elseif (FORCING)
    message("${Red} CFE BUILD 'FORCING' CASE!${ColourReset}")
    set(exe_name "cfe_forcing")
elseif (FORCINGPET)
    message("${Red} CFE BUILD 'FORCING AND PET' CASE!${ColourReset}")
    set(exe_name "cfe_forcingpet")
elseif (AETROOTZONE)
    message("${Red} CFE BUILD 'FORCING AND PET AND ROOTZONE-BASED AET' CASE!${ColourReset}")
    set(exe_name "cfe_aet_rootzone")
elseif (NGEN)
    message("${Red} CFE BUILD NGEN CASE!${ColourReset}")
elseif (UNITTEST)
    message("${Red} CFE BUILD UNITTEST CASE!${ColourReset}")
    set(exe_name "cfe_unittest")
endif ()

# set the project name
project(cfebmi VERSION 1.0.0 DESCRIPTION "OWP CFE BMI Module Shared Library")

IF (CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF (CMAKE_BUILD_TYPE MATCHES Debug)

message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# add the executable

## cfe + aorc + pet + smp
if (AETROOTZONE)
    add_executable(${exe_name} ./src/main_cfe_aorc_pet_rz_aet.cxx ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c
            ./src/conceptual_reservoir.c ./src/nash_cascade.c ./extern/aorc_bmi/src/aorc.c
            ./extern/aorc_bmi/src/bmi_aorc.c ./extern/evapotranspiration/src/pet.c
            ./extern/evapotranspiration/src/bmi_pet.c)

    add_library(cfelib ./extern/SoilMoistureProfiles/src/bmi_soil_moisture_profile.cxx
            ./extern/SoilMoistureProfiles/src/soil_moisture_profile.cxx
            ./extern/SoilMoistureProfiles/include/bmi_soil_moisture_profile.hxx
            ./extern/SoilMoistureProfiles/include/soil_moisture_profile.hxx)
    target_link_libraries(${exe_name} LINK_PUBLIC cfelib)
elseif (FORCING)
    add_executable(${exe_name} ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
            ./src/nash_cascade.c ./extern/aorc_bmi/src/aorc.c ./extern/aorc_bmi/src/bmi_aorc.c)
elseif (FORCINGPET)
    add_executable(${exe_name} ./src/main_cfe_aorc_pet.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
            ./src/nash_cascade.c ./extern/aorc_bmi/src/aorc.c ./extern/aorc_bmi/src/bmi_aorc.c
            ./extern/evapotranspiration/src/pet.c ./extern/evapotranspiration/src/bmi_pet.c)
elseif (BASE)
    add_executable(${exe_name} ./src/main.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
            ./src/nash_cascade.c)
elseif (UNITTEST)
    add_executable(${exe_name} test/combined_bmi_funcs_test.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
            ./src/nash_cascade.c)
endif ()


if (NOT NGEN)
    target_link_libraries(${exe_name} PRIVATE m)
    target_include_directories(${exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endif ()


# ngen CFE (add shared library)
if (NGEN)
    set(CFE_LIB_NAME_CMAKE cfebmi)
    set(CFE_LIB_DESC_CMAKE "OWP CFE BMI Module Shared Library")

    # Make sure these are compiled with this directive
    add_compile_definitions(BMI_ACTIVE)

    if (WIN32)
        add_library(cfebmi ./src/bmi_cfe.c ./src/cfe.c ./src/giuh.c ./src/conceptual_reservoir.c ./src/nash_cascade.c)
    else ()
        add_library(cfebmi SHARED ./src/bmi_cfe.c ./src/cfe.c ./src/giuh.c ./src/conceptual_reservoir.c ./src/nash_cascade.c)
    endif ()

    target_include_directories(cfebmi PRIVATE include)

    set_target_properties(cfebmi PROPERTIES VERSION ${PROJECT_VERSION})

    set_target_properties(cfebmi PROPERTIES PUBLIC_HEADER ./include/bmi_cfe.h)

    # Add unit testing using CMake and CTests, and custom test code
    enable_testing()

    ####################################################################################################################
    ####################################################################################################################
    # Include executable and (single) test for the legacy batch test file
    add_executable(combined_bmi_funcs_test test/combined_bmi_funcs_test.c)
    target_include_directories(combined_bmi_funcs_test PUBLIC ${PROJECT_SOURCE_DIR}/include test/include)
    target_link_libraries(combined_bmi_funcs_test cfebmi)
    add_test(NAME combined_bmi_funcs_test_1 COMMAND combined_bmi_funcs_test ${PROJECT_SOURCE_DIR}/configs/cfe_config_cat_87_pass.txt)
    # End combined_bmi_funcs_test tests
    ####################################################################################################################

    ####################################################################################################################
    ####################################################################################################################
    # Include executable and individual, isolated unit tests for BMI related functions
    # Control which test is run via custom implemented arg parsing within the executable
    add_executable(test_bmi_model ./test/test_bmi_model.c test/bmi_test_utils.c test/bmi_test_utils.c test/general_test_utils.c)
    target_include_directories(test_bmi_model PUBLIC ${PROJECT_SOURCE_DIR}/include test/include)
    target_link_libraries(test_bmi_model cfebmi)
    target_compile_definitions(test_bmi_model PUBLIC BMI_INIT_CONFIG_EX_1="${PROJECT_SOURCE_DIR}/configs/cfe_config_cat_87_pass.txt")
    #set(BMI_CONFIG_1 ${PROJECT_SOURCE_DIR}/configs/cfe_config_cat_87_pass.txt)

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
    # End test_bmi_model tests
    ####################################################################################################################

    include(GNUInstallDirs)

    install(TARGETS cfebmi
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    configure_file(cfebmi.pc.in cfebmi.pc @ONLY)

    install(FILES ${CMAKE_BINARY_DIR}/cfebmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

endif ()

unset(BASE CACHE)
unset(FORCING CACHE)
unset(FORCINGPET CACHE)
unset(AETROOTZONE CACHE)
unset(NGEN CACHE)
