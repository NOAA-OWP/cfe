# Build and Run Instructions
Detailed instructions on how to build and run CFE in three modes (standalone, pseudo, and nextgen frameworks) are provided below. Building CFE requires [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine.

## Clone the repository
```
git clone https://github.com/NOAA-OWP/cfe
cd cfe
git submodule update --init
mkdir build && cd build
```

**Notes:**
 - Before running the following examples, it is recommended to run the unittests [tests](https://github.com/NOAA-OWP/cfe/tree/master/test).
 - All build commands are run within the `build` directory, and run commands are executed from the cfe directory.
 
## Example 1 (standalone mode)
CFE reads local forcing data (standalone CFE BMI run), potential ET is not included. The example uses CAMELS catchment-87 data. **Notation:** BASE
### Build
```
cmake ../ -DBASE=ON
make && cd ..
```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/cfe/blob/master/run_cfe.sh">./run_cfe.sh</a> BASE   
</pre>

## Example 2 (pseudo framework mode)
This example couples CFE with AORC. AORC reads the forcing data and passes rain precipitation to CFE using BMI, potential ET is not included. The example uses CAMELS catchment-87 data. **Notation:** FORCING
### Build
```
cmake ../ -DFORCING=ON
make && cd ..
```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/cfe/blob/master/run_cfe.sh">./run_cfe.sh</a> FORCING
</pre>

## Example 3 (pseudo framework mode)
This example couples CFE with AORC and PET. AORC reads the forcing data and passes data to CFE and PET using BMIs. The actual ET is extracted directly from precipitation and from the soil using the Budyko function. The example uses CAMELS catchment-87 data. **Notation:** FORCINGPET
### Build
```
cmake ../ -DFORCINGPET=ON
make && cd ..
```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/cfe/blob/master/run_cfe.sh">./run_cfe.sh</a> FORCINGPET
</pre>

## Example 4 (pseudo framework mode)
The setup of this example is identical to example #3, however, actual ET is computed using a rootzone-based scheme. The example requires CFE coupling with [SoilMoistureProfiles](https://github.com/NOAA-OWP/SoilMoistureProfiles) module which provides a one-dimensional soil moisture profile that is used to estimate actual evapotranspiration (AET) from the deepest rootzone layer. **Notation:** AETROOTZONE

**Clone SoilMoistureProfiles repo:**
```
git clone https://github.com/NOAA-OWP/SoilMoistureProfiles extern/SoilMoistureProfiles (from cfe directory)
```
### Build
```
cmake ../ -DAETROOTZONE=ON
make && cd ..
```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/cfe/blob/master/run_cfe.sh">./run_cfe.sh</a> AETROOTZONE
</pre>

## Example 5 (nextgen framework mode)
This example couples CFE with PET. Nextgen framework reads the forcing data and passes data to CFE and PET using BMIs. The example uses CAMELS catchment-87 data. (Notation: NGEN)
### Build
See general instructions [here](https://github.com/NOAA-OWP/ngen/wiki/NGen-Tutorial) or [here](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/INSTALL.md#build-2) to build and run models in the ngen framework.
- ### Specific instructions for building nextgen and models needed for the example here
 - git clone https://github.com/noaa-owp/ngen && cd ngen
  - git submodule update --init --recursive
  - #### fortran bmi
    - cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
    - make -C extern/iso_c_fortran_bmi/cmake_build
  - #### build ngen
     - cmake -B cmake_build -S . -DBMI_C_LIB_ACTIVE=ON -DBMI_FORTRAN_ACTIVE=ON -DNGEN_ACTIVATE_PYTHON=ON
     - make -j4 -C cmake_build
  - #### CFE
    - git submodule update --remote extern/cfe/cfe 
    - cmake -B extern/cfe/cmake_build -S extern/cfe/cfe/ -DNGEN=ON
    - make -C extern/cfe/cmake_build
  - #### SLoTH
    SLoTH is also needed to run SFT in the ngen framework. SLoTH is a BMI that is used to set a bmi variable(s) that is not provided by other BMIs but required by the model. So build [SLoTH](https://github.com/NOAA-OWP/SLoTH) using the following instructions
    - cd extern/sloth/ && git checkout latest 
    - git submodule update --init --recursive
    - cd ../..
    - cmake -B extern/sloth/cmake_build -S extern/sloth/
    - make -C extern/sloth/cmake_build
### Run
The following pre-process step needs to be completed before running the example.
  #### Pre-process step
  ```
  mkdir cfe && cd cfe
  ln -s ../extern
  ln -s ../data
  ln -s ./extern/cfe/cfe/realizations
  ```
  
  **Note:** Make sure the "library_file" and "init_config" in the BMI blocks in the realization file are pointing to the right files, these paths depend on how you build your models.

  ```
  Run: ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_cfe_pet_surfgiuh.json
  ```

## Alternative: Compiling and Running CFE
### Example 1. Read local forcing file
To compile and run CFE with locally read forcing data, run the following from the command line:

1. `gcc -lm -Iinclude ./src/main.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c -o cfe_base`. This will generate an executable called `cfe_base`.
2.  Then run the model with example forcing data: `./cfe_base ./configs/cat_58_bmi_config_cfe.txt`  


### Example 2. CFE Model gets forcings passed from BMI

CFE was designed to read its own forcing file, but we have added an option to get forcings passed in through BMI using its `set_value` functionality. To demonstrate this functionality we have included the BMI-enabled AORC forcing read module. Follow the steps below:  

1. `gcc -lm -Iinclude ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c . ./extern/aorc_bmi/src/aorc.c ./extern/aorc_bmi/src/bmi_aorc.c  -o cfe_forcing`. This generates an executable called `cfe_forcing`. 

2. To run this executable you must pass the path to the corresponding configuration files for **BOTH** CFE and AORC (in that order): `./cfe_forcing ./configs/cat_89_bmi_config_cfe_pass.txt ./configs/cat_89_bmi_config_aorc.txt`

### Example 3. CFE Model gets forcings AND potential evapotranspiration passed from BMI
CFE can remove mass from the modeled system through evapotranspiration (directly from precipitation and from the soil using the Budyko function). Follow the steps below:  

1. `gcc -lm -Iinclude ./src/main_cfe_aorc_pet.c ./extern/evapotranspiration/src/pet.c ./extern/evapotranspiration/src/bmi_pet.c ./src/cfe.c ./src/bmi_cfe.c ./src/giuh.c ./src/conceptual_reservoir.c ./extern/aorc_bmi/src/aorc.c ./extern/aorc_bmi/src/bmi_aorc.c -o cfe_forcingpet`. This generates an executable called `cfe_forcingpet`.

2. To run this executable you must pass the path to the corresponding configuration files for CFE, PET and AORC (in that order):  `./cfe_forcingpet ./configs/cat_89_bmi_config_cfe_pass.txt ./configs/cat_89_bmi_config_aorc.txt ./configs/cat_89_bmi_config_pet_pass.txt`

### Example 4.
CFE rootzone-based example couples C and C++ modules and should be built with cmake, follow the instructions [here](https://github.com/NOAA-OWP/cfe/blob/master/INSTALL.md#example-4-pseudo-framework-mode).

### NOTES:
 - The configuration files must be passed in this order: (1) the CFE configuration file, (2) the forcing configuration file, (3) the potential evapotranspiration (PET) configuration file, and (4) the soil moisture profile configuration file
 - Original author code. The code includes a full program to read and process atmospheric forcing data, print the model output and check for mass balance closure. This code can be run from the [original_author_code](./original_author_code) directory. This code does not have a BMI implementation.
