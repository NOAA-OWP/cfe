import time
import numpy as np
import pandas as pd
import json
import cfe_orig
cfe1 = cfe_orig.CFE('./cat_58_config_cfe.json')
cfe1.run_unit_test()
