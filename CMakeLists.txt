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

if( (NOT BASE) AND (NOT FORCING) AND (NOT FORCINGPET) AND (NOT AETROOTZONE) AND (NOT NGEN) AND (NOT UNITTEST) )
  message("${Red}Options: BASE, FORCING, FORCINGPET, AETROOTZONE, NGEN, UNITTEST" ${ColourReset})
  message(FATAL_ERROR "Invalid option is provided, CMake will exit." )
endif()

if(BASE)
message("${Red} CFE BUILD 'BASE' CASE!${ColourReset}")
set(exe_name "cfe_base")
elseif(FORCING)
message("${Red} CFE BUILD 'FORCING' CASE!${ColourReset}")
set(exe_name "cfe_forcing")
elseif(FORCINGPET)
message("${Red} CFE BUILD 'FORCING AND PET' CASE!${ColourReset}")
set(exe_name "cfe_forcingpet")
elseif(AETROOTZONE)
message("${Red} CFE BUILD 'FORCING AND PET AND ROOTZONE-BASED AET' CASE!${ColourReset}")
set(exe_name "cfe_aet_rootzone")
elseif(NGEN)
message("${Red} CFE BUILD NGEN CASE!${ColourReset}")
elseif(UNITTEST)
message("${Red} CFE BUILD UNITTEST CASE!${ColourReset}")
set(exe_name "cfe_unittest")
endif()

# set the project name
project(cfebmi VERSION 1.0.0 DESCRIPTION "OWP CFE BMI Module Shared Library")

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# add the executable

## cfe + aorc + pet + smp
if(AETROOTZONE)
add_executable(${exe_name} ./src/main_cfe_aorc_pet_rz_aet.cxx ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c
			   ./src/conceptual_reservoir.c ./src/nash_cascade.c ./extern/aorc_bmi/src/aorc.c
			   ./extern/aorc_bmi/src/bmi_aorc.c ./extern/evapotranspiration/src/pet.c
			   ./extern/evapotranspiration/src/bmi_pet.c)

add_library(cfelib ./extern/SoilMoistureProfiles/src/bmi_soil_moisture_profile.cxx
		   ./extern/SoilMoistureProfiles/src/soil_moisture_profile.cxx
                   ./extern/SoilMoistureProfiles/include/bmi_soil_moisture_profile.hxx
		   ./extern/SoilMoistureProfiles/include/soil_moisture_profile.hxx)
target_link_libraries(${exe_name} LINK_PUBLIC cfelib)
elseif(FORCING)
add_executable(${exe_name} ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
			   ./src/nash_cascade.c ./extern/aorc_bmi/src/aorc.c ./extern/aorc_bmi/src/bmi_aorc.c)
elseif(FORCINGPET)
add_executable(${exe_name} ./src/main_cfe_aorc_pet.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
			   ./src/nash_cascade.c ./extern/aorc_bmi/src/aorc.c ./extern/aorc_bmi/src/bmi_aorc.c
			   ./extern/evapotranspiration/src/pet.c ./extern/evapotranspiration/src/bmi_pet.c)
elseif(BASE)
add_executable(${exe_name} ./src/main.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
			   ./src/nash_cascade.c)
elseif(UNITTEST)
add_executable(${exe_name} ./test/main_unit_test.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c
			   ./src/nash_cascade.c)
endif()


if(NOT NGEN)
target_link_libraries(${exe_name} PRIVATE m)
target_include_directories(${exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endif()


# ngen CFE (add shared library)
if(NGEN)
set(CFE_LIB_NAME_CMAKE cfebmi)
set(CFE_LIB_DESC_CMAKE "OWP CFE BMI Module Shared Library")

# Make sure these are compiled with this directive
add_compile_definitions(BMI_ACTIVE)

if(WIN32)
    add_library(cfebmi ./src/bmi_cfe.c ./src/cfe.c ./src/giuh.c ./src/conceptual_reservoir.c ./src/nash_cascade.c)
else()
    add_library(cfebmi SHARED ./src/bmi_cfe.c ./src/cfe.c ./src/giuh.c ./src/conceptual_reservoir.c ./src/nash_cascade.c)
endif()

target_include_directories(cfebmi PRIVATE include)

set_target_properties(cfebmi PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(cfebmi PROPERTIES PUBLIC_HEADER ./include/bmi_cfe.hxx)

include(GNUInstallDirs)

install(TARGETS cfebmi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(cfebmi.pc.in cfebmi.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/cfebmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

endif()

unset(BASE CACHE)
unset(FORCING CACHE)
unset(FORCINGPET CACHE)	
unset(AETROOTZONE CACHE)
unset(NGEN CACHE)
