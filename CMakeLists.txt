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

## cfe + aorc + pet + ftm
add_executable(cfe_smp ./src/main_cfe_aorc_pet_ftm.cxx ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c ./forcing_code/src/pet.c ./forcing_code/src/bmi_pet.c)

#add_executable(cfeft ./src/main_cfe_aorc_pet_ftm.cxx ./cfe/src/cfe.c ./cfe/src/bmi_cfe.c ./cfe/forcing_code/src/aorc.c ./cfe/forcing_code/src/bmi_aorc.c ./cfe/forcing_code/src/pet.c ./cfe/forcing_code/src/bmi_pet.c)

#add_executable(cfeft ./src/main_cfe_aorc_pet_ftm.cxx ./../cfe_owp/src/cfe.c ./../cfe_owp/src/bmi_cfe.c ./../cfe_owp/forcing_code/src/aorc.c ./../cfe_owp/forcing_code/src/bmi_aorc.c ./../cfe_owp/forcing_code/src/pet.c ./../cfe_owp/forcing_code/src/bmi_pet.c)

add_library(cfelib ./smc_coupler/src/bmi_coupler.cxx ./smc_coupler/src/smc_profile.cxx ./smc_coupler/include/bmi_coupler.hxx ./smc_coupler/include/smc_profile.hxx)

#add_library(freezethaw SHARED ./src/bmi_freezethaw.cxx ./src/freezethaw.cxx ./include/bmi_freezethaw.hxx ./include/freezethaw.hxx)
#add_library(cfebmi SHARED src/bmi_cfe.c src/cfe.c)

target_link_libraries(cfe_smp LINK_PUBLIC cfelib)
target_link_libraries(cfe_smp PRIVATE m)


##cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_BUILD_TYPE=Debug ../
##cmake --build .



