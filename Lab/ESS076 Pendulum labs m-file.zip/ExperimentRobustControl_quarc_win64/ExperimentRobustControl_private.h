/*
 * ExperimentRobustControl_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "ExperimentRobustControl".
 *
 * Model version              : 1.728
 * Simulink Coder version : 8.11 (R2016b) 25-Aug-2016
 * C source code generated on : Tue Aug 28 15:42:25 2018
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_ExperimentRobustControl_private_h_
#define RTW_HEADER_ExperimentRobustControl_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#include "ExperimentRobustControl.h"

/* A global buffer for storing error messages (defined in quanser_common library) */
EXTERN char _rt_error_message[512];
extern void Experime_PendulumAngleWatchdogX(RT_MODEL_ExperimentRobustCont_T *
  const ExperimentRobustControl_M, real_T rtu_alpharad,
  B_PendulumAngleWatchdogX_Expe_T *localB, P_PendulumAngleWatchdogX_Expe_T
  *localP);
extern void ExperimentR_SRV02XAngleWatchdog(RT_MODEL_ExperimentRobustCont_T *
  const ExperimentRobustControl_M, real_T rtu_thetarad,
  B_SRV02XAngleWatchdog_Experim_T *localB, P_SRV02XAngleWatchdog_Experim_T
  *localP);
extern void Experime_PendulumAngleWatchdogY(RT_MODEL_ExperimentRobustCont_T *
  const ExperimentRobustControl_M, real_T rtu_alphadeg,
  B_PendulumAngleWatchdogY_Expe_T *localB, P_PendulumAngleWatchdogY_Expe_T
  *localP);
extern void ExperimentR_SRV02YAngleWatchdog(RT_MODEL_ExperimentRobustCont_T *
  const ExperimentRobustControl_M, real_T rtu_thetarad,
  B_SRV02YAngleWatchdog_Experim_T *localB, P_SRV02YAngleWatchdog_Experim_T
  *localP);

/* private model entry point functions */
extern void ExperimentRobustControl_derivatives(void);

#endif                                 /* RTW_HEADER_ExperimentRobustControl_private_h_ */
