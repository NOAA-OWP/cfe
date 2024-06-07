#!/bin/bash
flag=$1

if [ ! $# == 1 ]; then
    echo "Usage: $0 OPTION={BASE,FORCING,FORCINGPET,AETROOTZONE}"
    echo "One of these options must be specified to run CFE"
  exit
fi

if [ $flag == "BASE" ] || [ "$flag" == "FORCING" ] || [ "$flag" == "FORCINGPET" ] || [ "$flag" == "AETROOTZONE" ]; then
    echo "CFE running with option $flag"
else
    echo "Invalid option! $flag"
    exit
fi


args=" "
exe_name=" "
if [ $flag == "BASE" ]; then
    args='./configs/cfe_config_cat_87.txt'
    exe_name='cfe_base'
else if [ $flag == "FORCING" ]; then
	 args="./configs/cfe_config_cat_87_pass.txt ./extern/aorc_bmi/configs/aorc_config_cat_87.txt"
	 exe_name='cfe_forcing'
     else if [ $flag == "FORCINGPET" ]; then
	      args='./configs/cfe_config_cat_87_pass.txt ./extern/aorc_bmi/configs/aorc_config_cat_87.txt ./extern/evapotranspiration/configs/pet_config_cat_87_pass.txt'
	      exe_name='cfe_forcingpet'
	  else if [ $flag == "AETROOTZONE" ]; then
		   args='./configs/cfe_config_laramie_pass_aet_rz.txt ./extern/aorc_bmi/configs/aorc_config_laramie.txt ./extern/evapotranspiration/configs/pet_config_laramie_pass.txt ./configs/smp_config_laramie.txt'
		   exe_name='cfe_aet_rootzone'
	       fi
	  fi
     fi
fi
echo "config file: $args"
./build/${exe_name} $args
