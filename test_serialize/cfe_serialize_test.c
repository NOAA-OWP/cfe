#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>   // for access()
#include "../include/bmi.h"
#include "../include/bmi_cfe.h"
#include "../include/serialize_state.h"
#include "../include/cfe_serialize_test.h"

/*
This main program creates two instances of a BMI-enabled version of
the CFE model.  It's purpose is to test a new set of BMI functions that
are intended to support model serialization for models written in C.

For Model 1, bmi.initialize() is called, and bmi.update() is called
several times to advance the model state.  Then serialize_state() is
called to retrieve the entire model state.

For Model 2, bmi.initialize() may be called, and then the function
deserialize_to_state() is called to put it in the same state as Model 1.

For both Models, bmi.update() is then called several more times and
results are compared with the compare_states() function.

In C, we cannot infer the type or size of an array, so we need BMI
functions to get: name, type, size and ptr.
Recall that array name is a pointer to first element.

Arrays are stored "flattened", so the model will take care of reading
values into 1D, 2D, 3D arrays.  We should only need the total number
of array elements.

Later on, ser_file should be an argument to main(), not hard-coded.
*/

//------------------------------------------------------------------------
int print_some(void *ptr_list[]){

  //------------------------------------------------
  // NOTE! Typecast ptr first, then add offset,
  //       into array, then dereference the ptr.
  //       Careful with order of operations.
  //------------------------------------------------
  puts("Printing some selected variables ...");

  printf("ptr_list[42] = volout = %f\n", *(double *)ptr_list[42]);
  printf("ptr_list[43] = volin  = %f\n", *(double *)ptr_list[43]);
  printf("ptr_list[46] = vol_in_nash    = %f\n", *(double *)ptr_list[46]);
  printf("ptr_list[47] = vol_out_nash   = %f\n", *(double *)ptr_list[47]);
  printf("ptr_list[51] = num_time_steps = %d\n", *(int *)ptr_list[51]);
  printf("ptr_list[52] = current_time_step = %d\n", *(int *)ptr_list[52]);
  printf("ptr_list[53] = time_step_size = %d\n", *(int *)ptr_list[53]);
  //----------------------------------------------
  // The following 3 methods all work, but don't
  // use "&" in get_state_var_ptrs().
  //----------------------------------------------
  printf("ptr_list[55] = forcing_file = %s\n", (char *)ptr_list[55]);
  //---------------------------------------------------------------------
  // printf("ptr_list[55] = forcing_file  = %s", ptr_list[55]);
  //---------------------------------------------------------------------
//   char *p;
//     // p =(char *)ptr_list[55];
//     p = ptr_list[55];
//     printf("ptr_list[55] = forcing_file  = ");
//     while (*p != '\0'){
//         printf("%c", *p++); } 
  //---------------------------------------------------------------
  printf("ptr_list[71] = aorc.time    = %ld\n", *(long *)ptr_list[71]);
  printf("ptr_list[75] = nash_storage = %f\n", *(double *)ptr_list[75]);
  //---------------------------------------------------------------    
  puts("");    // newline is added

  return 0;
}

//------------------------------------------------------------------------
// int main(void)
int main(int argc, const char *argv[])
{
  // Check for config file arg
  const char *cfg_file;
  if (argc<=1){
      cfg_file = "./configs/cat_89_bmi_config_cfe.txt";
      printf("WARNING: Missing config file argument.\n");
      printf("         Using default:\n");
      printf("         %s\n", cfg_file);
  } else {
      cfg_file = argv[1];
  }
  if( access( cfg_file, F_OK ) != 0 ) {
      printf("ERROR: cfg_file not found.\n");
      puts( cfg_file);
      puts("");
      exit(1);
  }
  //--------------------------------------------------------------
  const char *ser_file = "./test_serialize/model_state.ser";  // make arg later
  int n_steps1  = 10; // n_steps for Model1 before serializing
  int n_steps2  = 50; // n_steps for models after deserializing
  int verbose   = 1;
  int print_obj = 1;  // Set to 1 to print values after deserializing
  int n_state_vars;
  int result;
  int test_getters = 1;

  //--------------------------------------------------------------   
  if (verbose){
      puts(""); 
      puts("Allocating memory for BMI CFE instances 1 & 2 ...");
  }
  Bmi *model1 = (Bmi *) malloc(sizeof(Bmi));
  Bmi *model2 = (Bmi *) malloc(sizeof(Bmi));

  //--------------------------------------------------------------
  if (verbose){ puts("Registering CFE models 1 & 2 ..."); }

  register_bmi_cfe(model1);
  register_bmi_cfe(model2);

  //--------------------------------------------------------------
  if (verbose){ puts("Initializing CFE models 1 & 2 ..."); }

  model1->initialize(model1, cfg_file);
  model2->initialize(model2, cfg_file);

  //-------------------------------------------------------------- 
  if (test_getters){

      //#########################################
      // See "num_nash_lf = 2" in bmi_cfe.c
      // padded to 3 causes intermittent
      // mismatch for the 3rd value.
      // Also see: Get_var_length in bmi_cfe.c.
      //#########################################
      int length, itemsize, nbytes, index, grid;
      int count_all, count_input, count_output;
      char *name = "infiltration_depth_m";  // index = 2
      // char name[] = "infiltration_depth_m";  // Doesn't work.
      char *role_all    = "all";
      char *role_input  = "input_from_bmi";
      char *role_output = "output_to_bmi";
      //-------------------------------      
      //char *type;
      //char *role;
      //char *location;
      //------------------------------- 
      char type[BMI_MAX_TYPE_NAME];
      char role[BMI_MAX_ROLE_NAME];
      char location[BMI_MAX_LOCATION_NAME];
      char units[BMI_MAX_UNITS_NAME];
      //------------------------------- 
      puts("###############################################");      
      puts("Testing bmi.get_model_var_count()...");
      model1->get_model_var_count(model1, role_all, &count_all);
      printf("  role  = %s\n", role_all);
      printf("  count = %d\n", count_all);
      puts("Testing bmi.get_model_var_count()...");      
      model1->get_model_var_count(model1, role_input, &count_input);
      printf("  role  = %s\n", role_input);
      printf("  count = %d\n", count_input);
      puts("Testing bmi.get_model_var_count()...");      
      model1->get_model_var_count(model1, role_output, &count_output);
      printf("  role  = %s\n", role_output);
      printf("  count = %d\n", count_output);
      //--------------------------------------------
      // Alloc mem for names array role_all
      char **names_test_all = NULL;
      names_test_all = (char**) malloc (sizeof(char *) * count_all);
        for (int i=0; i<count_all; i++){
      names_test_all[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
      }
      // Populate array via get_model_var_names(role_all)
      puts("Testing bmi.get_model_var_names()...");
      printf("  role = %s\n", role_all);
      model1->get_model_var_names(model1, role_all, names_test_all);
      for (int j=0; j<count_all; j++){
          printf("  names_test_all[%d] = %s\n", j, names_test_all[j]);
      }
      // Free mem
      for (int n=0; n<count_all; n++){
          free (names_test_all[n]);
      }
      free (names_test_all);
      //--------------------------------------------
      // Alloc mem for names array role_input
      char **names_test_input = NULL;
      names_test_input = (char**) malloc (sizeof(char *) * count_input);
        for (int i=0; i<count_input; i++){
      names_test_input[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
      }
      // Populate array via get_model_var_names(role_input)
      puts("Testing bmi.get_model_var_names()...");
      printf("  role = %s\n", role_input);
      model1->get_model_var_names(model1, role_input, names_test_input);
      for (int j=0; j<count_input; j++){
          printf("  names_test_input[%d] = %s\n", j, names_test_input[j]);
      }
      // Free mem
      for (int n=0; n<count_input; n++){
          free (names_test_input[n]);
      }
      free (names_test_input);
      //--------------------------------------------
      // Alloc mem for names array role_output
      char **names_test_output = NULL;
      names_test_output = (char**) malloc (sizeof(char *) * count_output);
        for (int i=0; i<count_output; i++){
      names_test_output[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
      }
      // Populate array via get_model_var_names(role_output)
      puts("Testing bmi.get_model_var_names()...");
      printf("  role = %s\n", role_output);
      model1->get_model_var_names(model1, role_output, names_test_output);
      for (int j=0; j<count_output; j++){
          printf("  names_test_output[%d] = %s\n", j, names_test_output[j]);
      }
      // Free mem
      for (int n=0; n<count_output; n++){
          free (names_test_output[n]);
      }
      free (names_test_output);
      //--------------------------------------------
      puts(""); 
      printf("Let name = %s\n", name);
      //-------------------------------------------- 
      puts("Testing bmi.get_var_index()...");
      model1->get_var_index(model1, name, &index);
      printf("  index = %d\n", index);
      //-------------------------------------------- 
      puts("Testing bmi.get_var_grid()...");
      model1->get_var_grid(model1, name, &grid);
      printf("  grid = %d\n", grid);      
      //--------------------------------------------      
      //-------------------------------------------- 
      puts("Testing bmi.get_var_type()...");
      model1->get_var_type(model1, name, type);
      printf("  type = %s\n", type);
      //-------------------------------------------- 
      puts("Testing bmi.get_var_units()...");
      model1->get_var_units(model1, name, units);
      printf("  units = %s\n", units); 
      //--------------------------------------------
      puts("Testing bmi.get_var_role()...");
      model1->get_var_role(model1, name, role);
      printf("  role = %s\n", role);
      //--------------------------------------------
      puts("Testing bmi.get_var_length()...");
      model1->get_var_length(model1, name, &length);
      printf("  length = %d\n", length); 
      //--------------------------------------------
      // This one uses bmi.get_var_type().
      puts("Testing bmi.get_var_itemsize()...");
      model1->get_var_itemsize(model1, name, &itemsize);
      printf("  itemsize = %d\n", itemsize);
      //--------------------------------------------
      // This one uses bmi.get_var_type().
      puts("Testing bmi.get_var_nbytes()...");
      model1->get_var_nbytes(model1, name, &nbytes);
      printf("  nbytes = %d\n", nbytes);
      //--------------------------------------------
      puts("Testing bmi.get_var_location()...");
      model1->get_var_location(model1, name, location);
      printf("  location = %s\n", location);
      //--------------------------------------------
      puts("###############################################"); 
      puts("");                     
  }
  
  //--------------------------------------------------------------
  if (verbose){
      puts("Updating CFE model 1 ...");
      printf("  n_steps1 = %i \n", n_steps1);
      puts("");
  }

  for (int i=1; i<=n_steps1; i++){
      model1->update(model1);
  }

  //--------------------------------------------------------------  
  if (verbose){
      puts("Calling get_state_var_count() on CFE model 1 ...");
      //puts("Calling BMI.get_state_var_count() on CFE model 1 ...");
  }

  result = get_state_var_count(model1, &n_state_vars);
  if (result == BMI_FAILURE){
      puts("ERROR in get_state_var_count(); returning.");
      return BMI_FAILURE;
  }
  //model1->get_state_var_count(model1, &n_state_vars);  

  //--------------------------------------------------------------  
  if (verbose){
      puts("Calling get_state_var_names() on CFE model 1 ...");
      //puts("Calling BMI.get_state_var_names() on CFE model 1 ...");
  }

  char **names = NULL;
  names = (char**) malloc (sizeof(char *) * n_state_vars);
  for (int i=0; i<n_state_vars; i++){
    names[i] = (char*) malloc (sizeof(char) * BMI_MAX_VAR_NAME);
  }  

  result = get_state_var_names(model1, names);
  if (result == BMI_FAILURE){
      puts("ERROR in get_state_var_names(); returning.");
      return BMI_FAILURE;
  }
  //model1->get_state_var_names(model1, names); 
  
  //--------------------------------------------------------------  
  if (verbose){
      puts("Calling get_state_var_ptrs() on CFE model 1 ...");
      //puts("Calling BMI.get_state_var_ptrs() on CFE model 1 ...");

  }
  //---------------------------------------------
  // For testing:  All 3 print "8" on my MacPro
  //---------------------------------------------
  //printf("Size of void*    = %lu\n", sizeof(void*));
  //printf("Size of int*     = %lu\n", sizeof(int*));
  //printf("Size of double*  = %lu\n", sizeof(double*));
  //printf("\n");

  //--------------------------------------------------------------- 
  // See: https://stackoverflow.com/questions/7798383/
  //      array-of-pointers-to-multiple-types-c/7799543
  //--------------------------------------------------------------- 
  void *ptr_list[ n_state_vars ];
  result = get_state_var_ptrs(model1, names, ptr_list);
  if (result == BMI_FAILURE){
      puts("ERROR in get_state_var_ptrs(); returning.");
      return BMI_FAILURE;
  }
  //model1->get_state_var_ptrs(model1, ptr_list);

  if (verbose){ print_some( ptr_list ); }

  //--------------------------------------------------------------  
  if (verbose){ puts("Calling serialize() on CFE model 1 ..."); }

  //------------------------------------------------
  // Serialize Model1 state and save to:  ser_file
  //------------------------------------------------
  serialize( model1, ser_file );

  //--------------------------------------------------------------
  if (verbose){
      puts("Calling deserialize_to_state() on CFE model 2 ...");
  }

  //-----------------------------------------------
  // Deserialize Model1 state saved in "ser_file"
  // and set it as the new state of Model2
  //-----------------------------------------------
  deserialize_to_state( ser_file, model2, print_obj );
  
  //--------------------------------------------------------------
  if (verbose){
      puts("Updating BMI CFE model 2 ...");
      printf("n_steps2 = %i \n", n_steps2);
      puts("");
  }

  for (int i=1; i<=n_steps2; i++){
      model2->update(model2);
  }
 
  //--------------------------------------------------------------
  if (verbose){
      puts("Updating BMI CFE model 1 ...");
      printf("n_steps2 = %i \n", n_steps2);
      puts("");
  }

  for (int i=1; i<=n_steps2; i++){
      model1->update(model1);
  }

  //--------------------------------------------------------------  
  if (verbose){
      puts("Calling get_state_var_ptrs() on CFE model 1 ...");
      //puts("Calling BMI.get_state_var_ptrs() on CFE model 1 ...");
  }

  result = get_state_var_ptrs(model1, names, ptr_list);
  if (result == BMI_FAILURE){
      puts("ERROR in get_state_var_ptrs(); returning.");
      return BMI_FAILURE;
  }
  // model1->get_state_var_ptrs(model1, ptr_list);

  if (verbose){ print_some( ptr_list ); }
  
  //--------------------------------------------------------------  
  if (verbose){
      puts("Calling get_state_var_ptrs() on CFE model 2 ...");
      //puts("Calling BMI.get_state_var_ptrs() on CFE model 2 ...");
  }

  result = get_state_var_ptrs(model2, names, ptr_list);
  if (result == BMI_FAILURE){
      puts("ERROR in get_state_var_ptrs(); returning.");
      return BMI_FAILURE;
  }
  // model2->get_state_var_ptrs(model2, ptr_list);

  if (verbose){ print_some( ptr_list ); }

  //--------------------------------------------------------------
  if (verbose){
      puts("Freeing memory...");
  }
    
  //--------------------------------   
  // Free memory for string arrays
  //--------------------------------
  for (int i=0; i<n_state_vars; i++){
      free (names[i]);
      //free (types[i]);
  }
  free (names);
  //free (types);

  //--------------------------------------------------------------
  if (verbose){ puts("Comparing CFE model 1 & 2 state vars ..."); }
  
  compare_states( model1, model2 );
  
  //--------------------------------------------------------------
  if (verbose){ puts("Finalizing BMI CFE models 1 & 2 ..."); }

  model1->finalize(model1);
  model2->finalize(model2);  

  if (verbose){ 
      puts("Finished with serialization test.\n");
  }
  return 0;
}

//------------------------------------------------------------------------





