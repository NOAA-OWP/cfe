# BMI Unit Testing
New BMI components introduced are categorized as follows,  
- Model control functions (4)
- Model information functions (5)
- Variable information functions (6)
- Time functions (5)
- Variable getter and setter functions (5)
- Model grid functions (16)

We will fully examine functionality of all applicable definitions.

To run the BMI component unit test, simply run `./make_and_run_bmi_unit_test.sh` within this [test](./make_and_run_bmi_unit_test.sh) directory. 

The script uses a catchment-89 configuration found [here](./configs/cat_89_bmi_config_cfe.txt).
Note that the actual testing loop is much smaller than the number of time steps or end time generated via configuration file.

Recall that BMI guides interoperability for model-coupling, where model components (i.e. inputs and outputs) are easily shared amongst each other.
When testing outside of a true framework, we consider the behavior of BMI function definitions, rather than any expected values they produce.
