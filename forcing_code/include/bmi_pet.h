#ifndef BMI_PET_H
#define BMI_PET_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "../../include/bmi.h"
#include "pet.h"

Bmi* register_bmi_pet(Bmi *model);

pet_model * new_bmi_pet();


int read_file_line_counts_pet(const char* file_name, 
                                 int* line_count, 
                                 int* max_line_length);
int read_init_config_pet(pet_model* model, const char* config_file);

#if defined(__cplusplus)
}
#endif

#endif
