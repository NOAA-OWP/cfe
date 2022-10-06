#!/bin/bash
flag=$1

if [ ! $# == 1 ]; then
    echo "Usage: $0 OPTION={BASE,FORCING,FORCINGPET,AETROOTZONE}"
    echo "One of these options must be specified to run CFE"
  exit
fi

if [ $flag == "BASE" ] || [ "$flag" == "FORCING" ] || [ "$flag" == "FORCINGPET" ] || [ "$flag" == "AETROOTZONE" ]; then
    #or ($flag != "FORCING") or ($flag != "FORCINGPET") or ($flag != "AETROOTZONE") ]; then
echo "CFE running with option $flag"
else
echo "Invalid option! $flag"
exit
fi


args=" "
exe_name=" "
if [ $flag == "BASE" ]; then
    args='./configs/cat_87_bmi_config_cfe.txt'
    exe_name='cfe_base'
else if [ $flag == "FORCING" ]; then
	 args="./configs/cat_87_bmi_config_cfe_pass.txt ./configs/cat_87_bmi_config_aorc.txt"
	 exe_name='cfe_forcing'
     else if [ $flag == "FORCINGPET" ]; then
	      args='./configs/cat_87_bmi_config_cfe_pass.txt ./configs/cat_87_bmi_config_aorc.txt ./configs/cat_87_bmi_config_pet_pass.txt'
	      exe_name='cfe_forcingpet'
	  else if [ $flag == "AETROOTZONE" ]; then
		   args='./configs/laramie_bmi_config_cfe_pass_aet_rz.txt ./configs/laramie_bmi_config_aorc.txt ./configs/laramie_bmi_config_pet_pass.txt ./configs/laramie_bmi_config_smc_coupler.txt'
		   exe_name='cfe_aet_rootzone'
	       fi
	  fi
     fi
fi
echo "config file: $args"
./build/${exe_name} $args
