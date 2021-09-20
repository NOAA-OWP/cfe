# BMI Unit Testing
New BMI components introduced are categorized as follows,  
- Model control functions (4)
- Model information functions (5)
- Variable information functions (6)
- Time functions (5)
- Variable getter and setter functions (5)
- Model grid functions (16)

To run the BMI component unit test, 
```
$ cd test
$ make_and_run_bmi_unit_test
```
Console output will show results from `BMI_SUCCESS` or simply print `BMI_FAILURE`. Note that the actual testing loop is much smaller the the number of timsteps or end time defined via configuration file. 
Recall that BMI guides interoperability for model-coupling, where model components (i.e. inputs and outputs) are easily shared amongst each other.
When testing outside of a true framework, we consider the behavior of BMI function definitions, rather than any expected values they produce.