# Conceptual Functional Equivalent (CFE) Model

This model is designed to be a simplified model of the National Water Model, which is functionally equivalent.  CFE (Conceptual Functional Equivalent) is a conceptual model written by Fred Ogden and was not originally BMI compatible.  To see the original author code, please refer to the `original_author_code` directory.  For more information on the hypotheses and ideas underpinning the CFE model, see the [T-shirt Approximation of the National Water Model versions 1.2, 2.0, and 2.1](https://github.com/NOAA-OWP/cfe/edit/AET_rootzone/README.md#t-shirt-approximation-of-the-national-water-model-versions-12-20-and-21) section of this document.  The remainder of this document discusses the BMI enabled and expanded CFE model. 

## Build and Run Instructions
Detailed instructions on how to build and run the CFE can be found in the [INSTALL](https://github.com/NOAA-OWP/cfe/blob/ajk/doc_update/INSTALL.md) guide.
 - Test examples highlights
   - Unittest (see [tests](https://github.com/NOAA-OWP/cfe/blob/ajk/doc_update/test/README.md))
   - Example 1 (standalone mode): CFE reads local forcing data
   - Example 2 (pseudo framework mode): CFE coupled to AORC (AORC provides forcing data through BMI)
   - Example 3 (pseudo framework mode): CFE coupled to AORC (provides forcing data through BMI) and PET (provides potential evapotranspiration via BMI)
   - Example 4 (pseudo framework mode): Example #3 repeated with rootzone-based actual evapotranspiration
   - Example 5 (nextgen framework mode): CFE coupled to PET module
   
## Model Configuration File
A detailed description of the parameters for model configuration is provided [here](https://github.com/NOAA-OWP/cfe/tree/ajk/doc_update/configs/README.md).

## Getting help
For questions, please contact XYZ, the main maintainer of the repository.

## Known issues or raise an issue
We are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or to report a bug, please open an issue.

## Getting involved
See general instructions to contribute to the model development ([instructions](https://github.com/NOAA-OWP/cfe/blob/ajk/doc_update/CONTRIBUTING.md)) or simply fork the repository and submit a pull request.
