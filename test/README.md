<!-- TODO: this needs to be updated for new testing -->

# Existing Tests
Currently there are two "sets" of tests, which are configured for use with `CTest` directly within the CMake config:
* Combined test function exercising (basically) all BMI functions (using the `combined_bmi_funcs_test` [test executable](#test-executables))
  * The previously existing tests, testing the entire module running a single long function
* Individual BMI function unit tests (using the `test_bmi_model` [test executable](#test-executables))
  * A set of tests each focused on a single BMI function

# Running Tests

## Building
First, create or reload a CMake build directory initialized for ngen (i.e., with `-DNGEN=ON`).   You may have already done this, though if there have been changes to [CMakeLists.txt](../CMakeLists.txt) since then, you will need to run it again. E.g.:

```shell
# Assuming ./cfe is the local source directory for the repo, and ./cfe/build is the name of the CMake dir
cd cfe
cmake -B build -DNGEN=ON -S .
```

> [!IMPORTANT]
> The created build directory after `-B` - `build` here - may/can be named something else.   Adjust example commands here as necessary for that difference in your setup.

From there, while you could also build the testing targets individually, the simplest thing is just to build everything in the project:

```shell
cmake --build build
```

## Executing 

> [!WARNING]
> These commands assume the tests have already been built as described above.

There are a several ways to execute the tests.  One is to use the CMake `test` target to run all the tests.

```shell
cmake --build build --target test
```

Alternatively, you could use the CTest executable directly, which comes with CMake.  Again, this runs all the tests.

```shell
ctest --test-dir build 
```

See [the documentation for CMake and CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html#testing-using-ctest) for options on controlling which tests are run.


### Test Details

The first (currently only) test case uses a catchment-89 configuration found [here](../configs/cat_89_bmi_config_cfe_unit_test.txt).

In particular, for the older, combined tests, note that the actual testing loop is much smaller than the number of time steps or end time generated via configuration file or otherwise.

Recall that BMI guides interoperability for model-coupling, where model components (i.e. inputs and outputs) are easily shared amongst each other.
When testing outside of a true framework, we consider the behavior of BMI function definitions, rather than any expected values they produce.

# Test Code

Tests are implemented within the different test code files described below.

Tests are configured within the main [CMakeLists.txt](../CMakeLists.txt) file for the project.  See [this page](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html) in the CMake documentation for more complete details on how CMake facilitates testing with CTest.

## Files

It is recommended files in the test directory follow these conventions for code organization.

### combined_bmi_funcs_test.c
A special case: this is an older test code file running (essentially) all the module's BMI functions within a single, long testing routine.  This contains useful testing but does not isolate functionality being tested, and so is not a true unit test.

### *_utils.[c,h]
Reusable, helper utility structs and functions that aid in testing should be declared/defined in filed named `*_utils.h`/`*_utils.c`.  For example, `bmi_test_utils.h` declares a struct for a test fixture and several helper functions with repeatedly needed operations related to BMI.  The `general_test_utils.h` file includes more general testing util functions, including those to compare expected and actual values of different types.

### test_*.c
A unit test source file, containing a `main` function and used as the command of one or more defined unit tests.  This should be used by some `add_executable` command in[CMakeLists.txt](../CMakeLists.txt) to create an executable target in the build config, which in turn is used as the `COMMAND` argument for one or more `add_test` commands.

For example, `test_bmi_model.c` is the source file for unit tests of all the CFE BMI functions.  Its `main` function parses arguments to determine which individual unit test function to call and which test example to use to setup the test fixture `struct` (see [above](#bmi_test_util)) that is passed to the unit test function.

## Configuring More Tests

Again note that [this CMake documentation](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html#how-does-cmake-facilitate-testing) has more complete details on how to add tests.  

The basic process is:

* Add new test code and/or new test example configs as needed
* If new test code included a new test executable with its own `main()`, add this in [CMakeLists.txt](../CMakeLists.txt) with the `add_executable` CMake command and otherwise configured as needed (e.g., using `target_include_directories`, `target_link_libraries`, etc.)
* Define the new test in [CMakeLists.txt](../CMakeLists.txt) using the CMake `add_test` command
* Reload the CMake project (i.e., rerun the `cmake -B build -DNGEN=ON -S .` command [discussed above](#building))

### Example

```CMake
add_executable(test_bmi_model ./test/test_bmi_model.c test/bmi_test_utils.c test/bmi_test_utils.c test/general_test_utils.c)
target_include_directories(test_bmi_model PUBLIC ${PROJECT_SOURCE_DIR}/include test/include)
target_link_libraries(test_bmi_model cfebmi)
target_compile_definitions(test_bmi_model PUBLIC BMI_INIT_CONFIG_EX_1="${PROJECT_SOURCE_DIR}/configs/cfe_config_cat_87_pass.txt")
...
add_test(NAME test_initialize_1_a COMMAND test_bmi_model test_initialize)
add_test(NAME test_update_1 COMMAND test_bmi_model test_update)
```

