cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red         "${Esc}[31m")
# module setup options

#option(BASE "BASE" ON)

option(FORCING "FORCING" OFF)
option(FORCINGPET "FORCINGPET" OFF)
option(AETROOTZONE "AETROOTZONE" OFF)

if(FORCING)
message("${Red} CFE RUNNING WITH FORCING!${ColourReset}")
set(exe_name "cfe_forcing")
elseif(FORCINGPET)
message("${Red} CFE RUNNING WITH FORCING AND PET!${ColourReset}")
set(exe_name "cfe_forcingpet")
elseif(AETROOTZONE)
message("${Red} CFE RUNNING WITH FORCING AND PET AND ROOTZONE BASED AET TURNED ON!${ColourReset}")
set(exe_name "cfe_aet_rootzone")
else()
message("${Red}CFE RUNNING WITH BASE CASE!${ColourReset}")
set(exe_name "cfe_base")
endif()

message("generated exe = ${exe_name}")
# set the project name
project(cfe)

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
if(AETROOTZONE)
add_executable(${exe_name} ./src/main_cfe_aorc_pet_rz_aet.cxx ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./extern/forcing_code/src/aorc.c ./extern/forcing_code/src/bmi_aorc.c ./extern/evapotranspiration/src/pet.c ./extern/evapotranspiration/src/bmi_pet.c)

add_library(cfelib ./smc_coupler/src/bmi_soil_moisture_profile.cxx ./smc_coupler/src/soil_moisture_profile.cxx ./smc_coupler/include/bmi_soil_moisture_profile.hxx ./smc_coupler/include/soil_moisture_profile.hxx)
elseif(FORCING)
add_executable(${exe_name} ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./extern/forcing_code/src/aorc.c ./extern/forcing_code/src/bmi_aorc.c)
elseif(FORCINGPET)
add_executable(${exe_name} ./src/main_cfe_aorc_pet.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./extern/forcing_code/src/aorc.c ./extern/forcing_code/src/bmi_aorc.c ./extern/evapotranspiration/src/pet.c ./extern/evapotranspiration/src/bmi_pet.c)
else()
add_executable(${exe_name} ./src/main.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c)

endif()
if(AETROOTZONE)
target_link_libraries(${exe_name} LINK_PUBLIC cfelib)
endif()
target_link_libraries(${exe_name} PRIVATE m)


unset(BASE CACHE)
unset(FORCING CACHE)
unset(FORCINGPET CACHE)	
unset(AETROOTZONE CACHE)
