/*
 * bmi_cfe.h — BMI interface for CFE v3
 *
 * Model state is CFE_Model_Context (defined in cfe_types.h),
 * stored in Bmi.data and accessed via the CONTEXT() macro in bmi_cfe.c.
 */

#ifndef CFE_BMI_CFE_H
#define CFE_BMI_CFE_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "bmi.h"
#include "cfe_types.h"

/* Fill a caller-provided Bmi struct with CFE function pointers (ngen convention) */
Bmi* register_bmi_cfe(Bmi *model);

/* Allocate a zeroed CFE_Model_Context (populated by Initialize) */
CFE_Model_Context* new_bmi_cfe(void);

/* Finalize + free a heap-allocated Bmi struct */
void delete_bmi_cfe(Bmi *model);

#if defined(__cplusplus)
}
#endif

#endif /* CFE_BMI_CFE_H */
