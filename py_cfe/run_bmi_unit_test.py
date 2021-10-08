"""Run BMI Unit Testing.
Author: jgarrett
Modified by Jonathan Frame
Date: 10/6/2021"""

import os
import sys
import numpy as np

from pathlib import Path
import bmi_cfe # This is the BMI that we will be running


# setup a "success counter" for number of passing and failing bmi functions
# keep track of function def fails (vs function call)
pass_count = 0
fail_count = 0
var_name_counter = 0
fail_list = []

def bmi_except(fstring):
    """Prints message and updates counter and list

    Parameters
    ----------
    fstring : str
        Name of failing BMI function 
    """
    
    global fail_count, fail_list, var_name_counter
    print("**BMI ERROR** in " + fstring)
    if (var_name_counter == 0):
        fail_count += 1
        fail_list.append(fstring)

# Define config path
cfg_file=Path('./cat_58_config_cfe.json')
bmi=bmi_cfe.BMI_CFE(cfg_file)

print("\nBEGIN BMI UNIT TEST\n*******************\n");


if os.path.exists(cfg_file):
    print(" configuration found: " + str(cfg_file))
else:
    print(" no configuration found, exiting...")
    sys.exit()

#-------------------------------------------------------------------
# initialize()
try: 
    bmi.initialize()
    print(" initializing...");
    pass_count += 1
except:
    bmi_except('initialize()')

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# BMI: Model Information Functions
#-------------------------------------------------------------------
#-------------------------------------------------------------------
print("\nMODEL INFORMATION\n*****************")

#-------------------------------------------------------------------
# get_component_name()
try:
    print (" component name: " + bmi.get_component_name())
    pass_count += 1
except:
    bmi_except('get_component_name()')

#-------------------------------------------------------------------
# get_input_item_count()
try:
    print (" input item count: " + str(bmi.get_input_item_count()))
    pass_count += 1
except:
    bmi_except('get_input_item_count()')

#-------------------------------------------------------------------
# get_output_item_count()
try:
    print (" output item count: " + str(bmi.get_output_item_count()))
    pass_count += 1
except:
    bmi_except('get_output_item_count()')

#-------------------------------------------------------------------
# get_input_var_names()
try:    
    # only print statement if names exist
    test_get_input_var_names = bmi.get_input_var_names()
    if len(test_get_input_var_names) > 0:
        print (" input var names: ")
        for var_in in test_get_input_var_names:
            print ("  " + var_in)
    pass_count += 1
except:
    bmi_except('get_input_var_names()')

#-------------------------------------------------------------------
# get_input_var_names()
try:    
    # only print statement if out var list not null
    test_get_output_var_names =  bmi.get_output_var_names()
    if len(test_get_output_var_names) > 0:
        print (" output var names: ")
        for var_out in test_get_output_var_names:
            print ("  " + var_out)
    pass_count += 1
except:
    bmi_except('get_output_item_count()')
    

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# BMI: Variable Information Functions
#-------------------------------------------------------------------
#-------------------------------------------------------------------
print("\nVARIABLE INFORMATION\n********************")

for var_name in (bmi.get_output_var_names() + bmi.get_input_var_names()):  
    print (" " + var_name + ":")

    #-------------------------------------------------------------------
    # get_var_units()
    try: 
        print ("  units: " + bmi.get_var_units(var_name))
        if var_name_counter == 0:
            pass_count += 1
    except:
        bmi_except('get_var_units()')
    
    #-------------------------------------------------------------------
    # get_var_itemsize()
    try:
        print ("  itemsize: " + str(bmi.get_var_itemsize(var_name)))
        if var_name_counter == 0:
            pass_count += 1
    except:
        bmi_except('get_var_itemsize()')

    #-------------------------------------------------------------------
    # get_var_type()
    try:
        print ("  type: " + str(bmi.get_var_type(var_name)))
        if var_name_counter == 0:
            pass_count += 1
    except:
        bmi_except('get_var_type()')

    #-------------------------------------------------------------------
    # get_var_nbytes()
    try:
        print ("  nbytes: " + str(bmi.get_var_nbytes(var_name)))
        if var_name_counter == 0:
            pass_count += 1
    except:
        bmi_except('get_var_nbytes()')

    #-------------------------------------------------------------------
    # get_var_grid
    try:
        print ("  grid id: " + str(bmi.get_var_grid(var_name)))
        if var_name_counter == 0:
            pass_count += 1
    except:
        bmi_except('get_var_grid()')

    #-------------------------------------------------------------------
    # get_var_location
    try:
        print ("  location: " + bmi.get_var_location(var_name))
        if var_name_counter == 0:
            pass_count += 1
    except:
        bmi_except('get_var_location()')

    var_name_counter += 1

# reset back to zero
var_name_counter = 0

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# BMI: Time Functions
#-------------------------------------------------------------------
#-------------------------------------------------------------------
print("\nTIME INFORMATION\n****************")

#-------------------------------------------------------------------
# get_start_time()
try:
    print (" start time: " + str(bmi.get_start_time()))
    pass_count += 1
except:
    bmi_except('get_start_time()')

#-------------------------------------------------------------------
# get_end_time()
try:
    print (" end time: " + str(bmi.get_end_time()))
    pass_count += 1
except:
    bmi_except('get_end_time()')

#-------------------------------------------------------------------
# get_current_time()
try:
    print (" current time: " + str(bmi.get_current_time()))
    pass_count += 1
except:
    bmi_except('get_current_time()')

#-------------------------------------------------------------------
# get_time_step()
try:
    print (" time step: " + str(bmi.get_time_step()))
    pass_count += 1
except:
    bmi_except('get_time_step()')

#-------------------------------------------------------------------
# get_time_units()
try:
    print (" time units: " + bmi.get_time_units())
    pass_count += 1
except:
    bmi_except('get_time_units()')


#-------------------------------------------------------------------
#-------------------------------------------------------------------
# BMI: Model Grid Functions
#-------------------------------------------------------------------
#-------------------------------------------------------------------
print("\nGRID INFORMATION\n****************")
grid_id = 0 # there is only 1
print (" grid id: " + str(grid_id))

#-------------------------------------------------------------------
# get_grid_rank()
try:
    print ("  rank: " + str(bmi.get_grid_rank(grid_id)))
    pass_count += 1
except:
    bmi_except('get_grid_rank()')

#-------------------------------------------------------------------
# get_grid_size()
try:    
    print ("  size: " + str(bmi.get_grid_size(grid_id)))
    pass_count += 1
except:
    bmi_except('get_grid_size()')

#-------------------------------------------------------------------
# get_grid_type()    
try:
    print ("  type: " + bmi.get_grid_type(grid_id))
    pass_count += 1
except:
    bmi_except('get_grid_type()')    


#-------------------------------------------------------------------
#-------------------------------------------------------------------
# BMI: Variable Getter and Setter Functions
#-------------------------------------------------------------------
#-------------------------------------------------------------------    
print ("\nGET AND SET VALUES\n******************")

for var_name in (bmi.get_output_var_names() + bmi.get_input_var_names()):     
    print (" " + var_name + ":" )

    #-------------------------------------------------------------------
    # set_value()
    try:
        this_set_value = -99.0
        bmi.set_value(var_name, this_set_value)
        print ("  set value: " + str(this_set_value))

        if var_name_counter == 0: 
            pass_count += 1
    except:
        bmi_except('set_value()')

    #-------------------------------------------------------------------
    # get_value()
    try:
        that_get_value = bmi.get_value(var_name)
        # check if set_value() passed then see if get/set values match
        if 'set_value()' not in fail_list:
            if that_get_value == this_set_value:
                print ("  get value: " + str(that_get_value) + " (values match)")
            else: 
                print ("  get value: " + str(that_get_value) + " (values DO NOT match)")
        else:
            print ("  get value: " + str(that_get_value))      
        if var_name_counter == 0: 
            pass_count += 1
    except:
        bmi_except('get_value()')

    #-------------------------------------------------------------------
    # get_value_ptr()
    try:
        that_get_value_ptr = bmi.get_value_ptr(var_name)
        # check if set_value() passed then see if get/set values match
        if 'set_value()' not in fail_list:
            if that_get_value_ptr == this_set_value:
                print ("  get value ptr: " + str(that_get_value) + " (values match)")
            else: 
                print ("  get value ptr: " + str(that_get_value) + " (values DO NOT match)")                
        else:
            print ("  get value ptr: " + str(that_get_value_ptr))
        if var_name_counter == 0: 
            pass_count += 1
    except:
        bmi_except('get_value_ptr()')

    #-------------------------------------------------------------------
    # set_value_at_indices()   
    try:
        this_set_value_at_indices = -11.0
        bmi.set_value_at_indices(var_name,[0], this_set_value_at_indices)
        #print ("  set value at indices: -9.0, and got value:", bmi.get_value(var_name))
        print ("  set value at indices: " + str(this_set_value_at_indices)) 
        if var_name_counter == 0: 
            pass_count += 1
    except:
        bmi_except('set_value_at_indices()')

    #-------------------------------------------------------------------
    # get_value_at_indices()    
    try: 
        dest0 = np.empty(bmi.get_grid_size(0), dtype=float)
        that_get_value_at_indices = bmi.get_value_at_indices(var_name, dest0, [0])
        # check if set_value_at_indices() passed then see if get/set values match
        if 'set_value_at_indices()' not in fail_list:
            if that_get_value_at_indices == this_set_value_at_indices:
                print ("  get value at indices: " + str(that_get_value_at_indices) + " (values match)")
            else: 
                print ("  get value at indices: " + str(that_get_value_at_indices) + " (values DO NOT match)")                
        # prob worth while to bounce against set_value() if get_value_at_indices() failed..        
        elif 'set_value()' not in fail_list:
            if that_get_value_at_indices == this_set_value:
                print ("  get value at indices: " + str(that_get_value_at_indices) + " (values match)")
            else: 
                print ("  get value at indices: " + str(that_get_value_at_indices) + " (values DO NOT match)")               
        else:
            # JMFrame NOTE: converting a list/array to a string probably won't work
            #print ("  get value at indices: " + str(bmi.get_value_at_indices(var_name, dest0, [0])))
            
            print ("  get value at indices: ", that_get_value_at_indices)
        
        if var_name_counter == 0: 
            pass_count += 1
    except: 
        bmi_except('get_value_at_indices()')

    var_name_counter += 1

# set back to zero
var_name_counter = 0

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# BMI: Control Functions
#-------------------------------------------------------------------
#-------------------------------------------------------------------   
print ("\nCONTROL FUNCTIONS\n*****************")    
    
#-------------------------------------------------------------------
# update()
try:
    bmi.update()
    # go ahead and print time to show iteration
    # wrap another try/except incase get_current_time() failed
    try: 
        print (" updating...        time " + str(bmi.get_current_time()));
    except: 
        print (" updating...");
    pass_count += 1
except:
    bmi_except('update()')

#-------------------------------------------------------------------
# update_until()
try:
    bmi.update_until(100)
    # go ahead and print time to show iteration
    # wrap another try/except incase get_current_time() failed
    try: 
        print (" updating until...  time " + str(bmi.get_current_time()))
    except: 
        print (" updating until...")
    pass_count += 1
except:
    bmi_except('update_until()')          

#-------------------------------------------------------------------
# finalize()
try:
    bmi.finalize()
    print (" finalizing...")
    pass_count += 1
except:
    bmi_except('finalize()')

# lastly - print test summary
print ("\n Total BMI function PASS: " + str(pass_count))
print (" Total BMI function FAIL: " + str(fail_count))
for ff in fail_list:
    print ("  " + ff)   