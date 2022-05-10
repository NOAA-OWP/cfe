cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})


# set the project name
project(cfe_smp)

set(CMAKE_BUILD_TYPE Debug)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# add the executable

## cfe + aorc + pet + smp
add_executable(cfe_smp ./src/main_cfe_aorc_pet_rzAET.cxx ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c ./forcing_code/src/pet.c ./forcing_code/src/bmi_pet.c)

add_library(cfelib ./smc_coupler/src/bmi_soil_moisture_profile.cxx ./smc_coupler/src/soil_moisture_profile.cxx ./smc_coupler/include/bmi_soil_moisture_profile.hxx ./smc_coupler/include/soil_moisture_profile.hxx)

target_link_libraries(cfe_smp LINK_PUBLIC cfelib)
target_link_libraries(cfe_smp PRIVATE m)




