#ifndef BMI_AORC_H
#define BMI_AORC_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "../../include/bmi.h"
#include "aorc.h"

Bmi* register_bmi_aorc(Bmi *model);

aorc_model * new_bmi_aorc();


int read_file_line_counts_aorc(const char* file_name, 
                                 int* line_count, 
                                 int* max_line_length);
int read_init_config_aorc(aorc_model* model, const char* config_file);

#if defined(__cplusplus)
}
#endif

#endif
