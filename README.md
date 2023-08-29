# Conceptual Functional Equivalent (CFE) Model

CFE (Conceptual Functional Equivalent) is a simplified conceptual model written by Fred Ogden that is designed to be functionally equivalent to the National Water Model. To see the original author code, which is not BMI compatible, please refer to the [original_author_code](https://github.com/NOAA-OWP/cfe/tree/master/original_author_code) directory.  For more information on the hypotheses and ideas underpinning the CFE model, see the [T-shirt Approximation of the National Water Model versions 1.2, 2.0, and 2.1](https://github.com/NOAA-OWP/cfe/blob/master/MODEL.md) section of this document.  The remainder of this document discusses the BMI enabled and expanded CFE model. 

## Build and Run Instructions
Detailed instructions on how to build and run CFE can be found in the [INSTALL](https://github.com/NOAA-OWP/cfe/blob/master/INSTALL.md) guide.
 - Test examples highlights
   - Unittest (see [tests](https://github.com/NOAA-OWP/cfe/blob/master/test/README.md))
   - Example 1 (standalone mode): CFE reads local forcing data
   - Example 2 (pseudo framework mode): CFE coupled to AORC (AORC provides forcing data through BMI)
   - Example 3 (pseudo framework mode): CFE coupled to AORC (provides forcing data through BMI) and PET (provides potential evapotranspiration via BMI)
   - Example 4 (pseudo framework mode): Example #3 repeated with rootzone-based actual evapotranspiration
   - Example 5 (nextgen framework mode): CFE coupled to PET module
   
## Model Configuration File
A detailed description of the parameters for model configuration is provided [here](https://github.com/NOAA-OWP/cfe/tree/master/configs/README.md).

## Getting help
For questions, please contact XYZ, the main maintainer of the repository.

## Known issues or raise an issue
We are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or to report a bug, please open an issue.

## Getting involved
See general instructions to contribute to the model development ([instructions](https://github.com/NOAA-OWP/cfe/blob/master/CONTRIBUTING.md)) or simply fork the repository and submit a pull request.
