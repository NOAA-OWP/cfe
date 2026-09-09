#ifndef CFE_BMI_SERIALIZATION
#define CFE_BMI_SERIALIZATION

#ifdef __cplusplus
extern "C" {
#endif

#include "bmi.h"

int free_serialized_cfe(Bmi* bmi);
int load_serialized_cfe(Bmi* bmi, char* data);
int new_serialized_cfe(Bmi* bmi);

#ifdef __cplusplus
}
#endif

#endif
