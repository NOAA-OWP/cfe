import time
import numpy as np
import pandas as pd
import json
import bmi_cfe
cfe1 = bmi_cfe.BMI_CFE('./cat_58_config_cfe.json')
cfe1.initialize()
cfe1.run_unit_test()
