/*
 * ExperimentRobustControl.c
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

#include "ExperimentRobustControl.h"
#include "ExperimentRobustControl_private.h"
#include "ExperimentRobustControl_dt.h"

/* Block signals (auto storage) */
B_ExperimentRobustControl_T ExperimentRobustControl_B;

/* Continuous states */
X_ExperimentRobustControl_T ExperimentRobustControl_X;

/* Block states (auto storage) */
DW_ExperimentRobustControl_T ExperimentRobustControl_DW;

/* Real-time model */
RT_MODEL_ExperimentRobustCont_T ExperimentRobustControl_M_;
RT_MODEL_ExperimentRobustCont_T *const ExperimentRobustControl_M =
  &ExperimentRobustControl_M_;

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 64;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  ExperimentRobustControl_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * Output and update for atomic system:
 *    '<S6>/Pendulum Angle Watchdog: X'
 *    '<S23>/Pendulum Angle Watchdog: X'
 */
void Experime_PendulumAngleWatchdogX(RT_MODEL_ExperimentRobustCont_T * const
  ExperimentRobustControl_M, real_T rtu_alpharad,
  B_PendulumAngleWatchdogX_Expe_T *localB, P_PendulumAngleWatchdogX_Expe_T
  *localP)
{
  /* Abs: '<S9>/Abs' */
  localB->Abs = fabs(rtu_alpharad);

  /* RelationalOperator: '<S9>/Check Maximum Allowed Angle' incorporates:
   *  Constant: '<S9>/X MAX'
   */
  localB->CheckMaximumAllowedAngle = (ExperimentRobustControl_P.ALPHA_MAX <
    localB->Abs);

  /* RelationalOperator: '<S12>/Compare' incorporates:
   *  Constant: '<S12>/Constant'
   */
  localB->Compare = (localB->CheckMaximumAllowedAngle != localP->Constant_Value);

  /* Stop: '<S11>/Stop Simulation' */
  if (localB->Compare) {
    rtmSetStopRequested(ExperimentRobustControl_M, 1);
  }

  /* End of Stop: '<S11>/Stop Simulation' */
}

/*
 * Output and update for atomic system:
 *    '<S6>/SRV02 X Angle Watchdog'
 *    '<S23>/SRV02 X Angle Watchdog'
 */
void ExperimentR_SRV02XAngleWatchdog(RT_MODEL_ExperimentRobustCont_T * const
  ExperimentRobustControl_M, real_T rtu_thetarad,
  B_SRV02XAngleWatchdog_Experim_T *localB, P_SRV02XAngleWatchdog_Experim_T
  *localP)
{
  /* Abs: '<S10>/Abs' */
  localB->Abs = fabs(rtu_thetarad);

  /* RelationalOperator: '<S10>/Check Maximum Allowed Angle' incorporates:
   *  Constant: '<S10>/X MAX'
   */
  localB->CheckMaximumAllowedAngle = (ExperimentRobustControl_P.THETA_MAX <
    localB->Abs);

  /* RelationalOperator: '<S14>/Compare' incorporates:
   *  Constant: '<S14>/Constant'
   */
  localB->Compare = (localB->CheckMaximumAllowedAngle != localP->Constant_Value);

  /* Stop: '<S13>/Stop Simulation' */
  if (localB->Compare) {
    rtmSetStopRequested(ExperimentRobustControl_M, 1);
  }

  /* End of Stop: '<S13>/Stop Simulation' */
}

/*
 * Output and update for atomic system:
 *    '<S7>/Pendulum Angle Watchdog: Y'
 *    '<S24>/Pendulum Angle Watchdog: Y'
 */
void Experime_PendulumAngleWatchdogY(RT_MODEL_ExperimentRobustCont_T * const
  ExperimentRobustControl_M, real_T rtu_alphadeg,
  B_PendulumAngleWatchdogY_Expe_T *localB, P_PendulumAngleWatchdogY_Expe_T
  *localP)
{
  /* Abs: '<S15>/Abs' */
  localB->Abs = fabs(rtu_alphadeg);

  /* RelationalOperator: '<S15>/Check Maximum Allowed Angle' incorporates:
   *  Constant: '<S15>/X MAX'
   */
  localB->CheckMaximumAllowedAngle = (ExperimentRobustControl_P.ALPHA_MAX <
    localB->Abs);

  /* RelationalOperator: '<S18>/Compare' incorporates:
   *  Constant: '<S18>/Constant'
   */
  localB->Compare = (localB->CheckMaximumAllowedAngle != localP->Constant_Value);

  /* Stop: '<S17>/Stop Simulation' */
  if (localB->Compare) {
    rtmSetStopRequested(ExperimentRobustControl_M, 1);
  }

  /* End of Stop: '<S17>/Stop Simulation' */
}

/*
 * Output and update for atomic system:
 *    '<S7>/SRV02 Y Angle Watchdog'
 *    '<S24>/SRV02 Y Angle Watchdog'
 */
void ExperimentR_SRV02YAngleWatchdog(RT_MODEL_ExperimentRobustCont_T * const
  ExperimentRobustControl_M, real_T rtu_thetarad,
  B_SRV02YAngleWatchdog_Experim_T *localB, P_SRV02YAngleWatchdog_Experim_T
  *localP)
{
  /* Abs: '<S16>/Abs' */
  localB->Abs = fabs(rtu_thetarad);

  /* RelationalOperator: '<S16>/Check Maximum Allowed Angle' incorporates:
   *  Constant: '<S16>/X MAX'
   */
  localB->CheckMaximumAllowedAngle = (ExperimentRobustControl_P.THETA_MAX <
    localB->Abs);

  /* RelationalOperator: '<S20>/Compare' incorporates:
   *  Constant: '<S20>/Constant'
   */
  localB->Compare = (localB->CheckMaximumAllowedAngle != localP->Constant_Value);

  /* Stop: '<S19>/Stop Simulation' */
  if (localB->Compare) {
    rtmSetStopRequested(ExperimentRobustControl_M, 1);
  }

  /* End of Stop: '<S19>/Stop Simulation' */
}

/* Model output function */
void ExperimentRobustControl_output(void)
{
  int_T ci;
  real_T y;
  real_T y_0[5];
  int32_T i;
  real_T u0;
  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
    /* set solver stop time */
    if (!(ExperimentRobustControl_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&ExperimentRobustControl_M->solverInfo,
                            ((ExperimentRobustControl_M->Timing.clockTickH0 + 1)
        * ExperimentRobustControl_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&ExperimentRobustControl_M->solverInfo,
                            ((ExperimentRobustControl_M->Timing.clockTick0 + 1) *
        ExperimentRobustControl_M->Timing.stepSize0 +
        ExperimentRobustControl_M->Timing.clockTickH0 *
        ExperimentRobustControl_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(ExperimentRobustControl_M)) {
    ExperimentRobustControl_M->Timing.t[0] = rtsiGetT
      (&ExperimentRobustControl_M->solverInfo);
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(ExperimentRobustControl_DW.ResetEncoders_SubsysRanBC);
  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S3>/HIL Read Encoder Timebase' */

    /* S-Function Block: ExperimentRobustControl/SRV02+2DIP-E/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (ExperimentRobustControl_DW.HILReadEncoderTimebase_Task, 1,
         &ExperimentRobustControl_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
      } else {
        ExperimentRobustControl_B.HILReadEncoderTimebase_o1 =
          ExperimentRobustControl_DW.HILReadEncoderTimebase_Buffer[0];
        ExperimentRobustControl_B.HILReadEncoderTimebase_o2 =
          ExperimentRobustControl_DW.HILReadEncoderTimebase_Buffer[1];
        ExperimentRobustControl_B.HILReadEncoderTimebase_o3 =
          ExperimentRobustControl_DW.HILReadEncoderTimebase_Buffer[2];
        ExperimentRobustControl_B.HILReadEncoderTimebase_o4 =
          ExperimentRobustControl_DW.HILReadEncoderTimebase_Buffer[3];
      }
    }

    /* Gain: '<S23>/2D IP Encoder Calibration (rad//count)' */
    y = ExperimentRobustControl_P.K_ENC_2DIP[0];
    ExperimentRobustControl_B.uDIPEncoderCalibrationradcount = y *
      ExperimentRobustControl_B.HILReadEncoderTimebase_o3;

    /* Outputs for Atomic SubSystem: '<S23>/Pendulum Angle Watchdog: X' */
    Experime_PendulumAngleWatchdogX(ExperimentRobustControl_M,
      ExperimentRobustControl_B.uDIPEncoderCalibrationradcount,
      &ExperimentRobustControl_B.PendulumAngleWatchdogX,
      (P_PendulumAngleWatchdogX_Expe_T *)
      &ExperimentRobustControl_P.PendulumAngleWatchdogX);

    /* End of Outputs for SubSystem: '<S23>/Pendulum Angle Watchdog: X' */

    /* Gain: '<S23>/Encoder Calibration  (rad//count)' */
    ExperimentRobustControl_B.EncoderCalibrationradcount =
      ExperimentRobustControl_P.K_ENC *
      ExperimentRobustControl_B.HILReadEncoderTimebase_o1;

    /* Outputs for Atomic SubSystem: '<S23>/SRV02 X Angle Watchdog' */
    ExperimentR_SRV02XAngleWatchdog(ExperimentRobustControl_M,
      ExperimentRobustControl_B.EncoderCalibrationradcount,
      &ExperimentRobustControl_B.SRV02XAngleWatchdog,
      (P_SRV02XAngleWatchdog_Experim_T *)
      &ExperimentRobustControl_P.SRV02XAngleWatchdog);

    /* End of Outputs for SubSystem: '<S23>/SRV02 X Angle Watchdog' */

    /* Gain: '<S24>/2D IP Encoder Y  Calibration (rad//count)' */
    y = ExperimentRobustControl_P.K_ENC_2DIP[1];
    ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount = y *
      ExperimentRobustControl_B.HILReadEncoderTimebase_o4;

    /* Outputs for Atomic SubSystem: '<S24>/Pendulum Angle Watchdog: Y' */
    Experime_PendulumAngleWatchdogY(ExperimentRobustControl_M,
      ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount,
      &ExperimentRobustControl_B.PendulumAngleWatchdogY,
      (P_PendulumAngleWatchdogY_Expe_T *)
      &ExperimentRobustControl_P.PendulumAngleWatchdogY);

    /* End of Outputs for SubSystem: '<S24>/Pendulum Angle Watchdog: Y' */

    /* Gain: '<S24>/Encoder Calibration Y (rad//count)' */
    ExperimentRobustControl_B.EncoderCalibrationYradcount =
      ExperimentRobustControl_P.K_ENC *
      ExperimentRobustControl_B.HILReadEncoderTimebase_o2;

    /* Outputs for Atomic SubSystem: '<S24>/SRV02 Y Angle Watchdog' */
    ExperimentR_SRV02YAngleWatchdog(ExperimentRobustControl_M,
      ExperimentRobustControl_B.EncoderCalibrationYradcount,
      &ExperimentRobustControl_B.SRV02YAngleWatchdog,
      (P_SRV02YAngleWatchdog_Experim_T *)
      &ExperimentRobustControl_P.SRV02YAngleWatchdog);

    /* End of Outputs for SubSystem: '<S24>/SRV02 Y Angle Watchdog' */
  }

  /* StateSpace: '<S1>/State-Space' */
  for (i = 0; i < 2; i++) {
    ExperimentRobustControl_B.StateSpace[i] = 0.0;
    for (ci = 0; ci < 54; ci++) {
      ExperimentRobustControl_B.StateSpace[i] += ExperimentRobustControl_P.ch
        [(ci << 1) + i] * ExperimentRobustControl_X.StateSpace_CSTATE[ci];
    }
  }

  /* End of StateSpace: '<S1>/State-Space' */

  /* Integrator: '<S6>/Integrator' */
  ExperimentRobustControl_B.Integrator_m =
    ExperimentRobustControl_X.Integrator_CSTATE_c;

  /* Integrator: '<S4>/Integrator' */
  ExperimentRobustControl_B.Integrator[0] =
    ExperimentRobustControl_X.Integrator_CSTATE[0];

  /* Integrator: '<S5>/Integrator' */
  ExperimentRobustControl_B.Integrator_k[0] =
    ExperimentRobustControl_X.Integrator_CSTATE_n[0];

  /* Integrator: '<S4>/Integrator' */
  ExperimentRobustControl_B.Integrator[1] =
    ExperimentRobustControl_X.Integrator_CSTATE[1];

  /* Integrator: '<S5>/Integrator' */
  ExperimentRobustControl_B.Integrator_k[1] =
    ExperimentRobustControl_X.Integrator_CSTATE_n[1];

  /* Integrator: '<S4>/Integrator' */
  ExperimentRobustControl_B.Integrator[2] =
    ExperimentRobustControl_X.Integrator_CSTATE[2];

  /* Integrator: '<S5>/Integrator' */
  ExperimentRobustControl_B.Integrator_k[2] =
    ExperimentRobustControl_X.Integrator_CSTATE_n[2];

  /* Integrator: '<S4>/Integrator' */
  ExperimentRobustControl_B.Integrator[3] =
    ExperimentRobustControl_X.Integrator_CSTATE[3];

  /* Integrator: '<S5>/Integrator' */
  ExperimentRobustControl_B.Integrator_k[3] =
    ExperimentRobustControl_X.Integrator_CSTATE_n[3];

  /* Integrator: '<S7>/Integrator' */
  ExperimentRobustControl_B.Integrator_c =
    ExperimentRobustControl_X.Integrator_CSTATE_a;

  /* SignalConversion: '<S2>/TmpSignal ConversionAtLQR Gain: YInport1' */
  ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[0] =
    ExperimentRobustControl_B.EncoderCalibrationYradcount;
  ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[1] =
    ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount;
  ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[2] =
    ExperimentRobustControl_B.Integrator_k[2];
  ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[3] =
    ExperimentRobustControl_B.Integrator_k[3];
  ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[4] =
    ExperimentRobustControl_B.Integrator_c;

  /* ManualSwitch: '<Root>/Manual Switch2' */
  if (ExperimentRobustControl_P.ManualSwitch2_CurrentSetting == 1) {
    ExperimentRobustControl_B.ManualSwitch2[0] =
      ExperimentRobustControl_B.StateSpace[0];
    ExperimentRobustControl_B.ManualSwitch2[1] =
      ExperimentRobustControl_B.StateSpace[1];
  } else {
    /* SignalConversion: '<S2>/TmpSignal ConversionAtLQR Gain: XInport1' */
    ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn[0] =
      ExperimentRobustControl_B.EncoderCalibrationradcount;
    ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn[1] =
      ExperimentRobustControl_B.uDIPEncoderCalibrationradcount;
    ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn[2] =
      ExperimentRobustControl_B.Integrator[2];
    ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn[3] =
      ExperimentRobustControl_B.Integrator[3];
    ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn[4] =
      ExperimentRobustControl_B.Integrator_m;

    /* Gain: '<S2>/LQR Gain: X' */
    u0 = 0.0;
    for (i = 0; i < 5; i++) {
      y = -ExperimentRobustControl_P.K_CT_labb[i];
      u0 += y * ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn[i];

      /* Gain: '<S2>/LQR Gain: Y' */
      y = -ExperimentRobustControl_P.K_CT_labb[i];
      y_0[i] = y;
    }

    ExperimentRobustControl_B.LQRGainX = u0;

    /* End of Gain: '<S2>/LQR Gain: X' */

    /* Gain: '<S2>/LQR Gain: Y' */
    u0 = 0.0;
    for (i = 0; i < 5; i++) {
      u0 += y_0[i] * ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[i];
    }

    ExperimentRobustControl_B.LQRGainY = u0;
    ExperimentRobustControl_B.ManualSwitch2[0] =
      ExperimentRobustControl_B.LQRGainX;
    ExperimentRobustControl_B.ManualSwitch2[1] =
      ExperimentRobustControl_B.LQRGainY;
  }

  /* End of ManualSwitch: '<Root>/Manual Switch2' */

  /* Saturate: '<S21>/Amplifier  Saturation (V)' */
  y = -ExperimentRobustControl_P.VMAX_AMP;

  /* Gain: '<Root>/Gain1' */
  ExperimentRobustControl_B.Gain1[0] = ExperimentRobustControl_P.OCL *
    ExperimentRobustControl_B.ManualSwitch2[0];

  /* Gain: '<S21>/Direction Convention: (Right-Hand) system' */
  ExperimentRobustControl_B.DirectionConventionRightHandsys[0] =
    ExperimentRobustControl_P.DirectionConventionRightHandsys *
    ExperimentRobustControl_B.Gain1[0];

  /* Saturate: '<S21>/Amplifier  Saturation (V)' */
  u0 = ExperimentRobustControl_B.DirectionConventionRightHandsys[0];
  if (u0 > ExperimentRobustControl_P.VMAX_AMP) {
    u0 = ExperimentRobustControl_P.VMAX_AMP;
  } else {
    if (u0 < y) {
      u0 = y;
    }
  }

  ExperimentRobustControl_B.AmplifierSaturationV[0] = u0;

  /* Gain: '<Root>/Gain1' */
  ExperimentRobustControl_B.Gain1[1] = ExperimentRobustControl_P.OCL *
    ExperimentRobustControl_B.ManualSwitch2[1];

  /* Gain: '<S21>/Direction Convention: (Right-Hand) system' */
  ExperimentRobustControl_B.DirectionConventionRightHandsys[1] =
    ExperimentRobustControl_P.DirectionConventionRightHandsys *
    ExperimentRobustControl_B.Gain1[1];

  /* Saturate: '<S21>/Amplifier  Saturation (V)' */
  u0 = ExperimentRobustControl_B.DirectionConventionRightHandsys[1];
  if (u0 > ExperimentRobustControl_P.VMAX_AMP) {
    u0 = ExperimentRobustControl_P.VMAX_AMP;
  } else {
    if (u0 < y) {
      u0 = y;
    }
  }

  ExperimentRobustControl_B.AmplifierSaturationV[1] = u0;

  /* Gain: '<S21>/Inverse Amplifier  Gain (V//V)' */
  y = ExperimentRobustControl_P.K_AMP;
  y = 1.0 / y;
  ExperimentRobustControl_B.InverseAmplifierGainVV[0] = y *
    ExperimentRobustControl_B.AmplifierSaturationV[0];
  ExperimentRobustControl_B.InverseAmplifierGainVV[1] = y *
    ExperimentRobustControl_B.AmplifierSaturationV[1];

  /* Saturate: '<S21>/DACB Saturation (V)' */
  y = -ExperimentRobustControl_P.VMAX_DAC;
  u0 = ExperimentRobustControl_B.InverseAmplifierGainVV[0];
  if (u0 > ExperimentRobustControl_P.VMAX_DAC) {
    u0 = ExperimentRobustControl_P.VMAX_DAC;
  } else {
    if (u0 < y) {
      u0 = y;
    }
  }

  ExperimentRobustControl_B.DACBSaturationV[0] = u0;
  u0 = ExperimentRobustControl_B.InverseAmplifierGainVV[1];
  if (u0 > ExperimentRobustControl_P.VMAX_DAC) {
    u0 = ExperimentRobustControl_P.VMAX_DAC;
  } else {
    if (u0 < y) {
      u0 = y;
    }
  }

  ExperimentRobustControl_B.DACBSaturationV[1] = u0;

  /* End of Saturate: '<S21>/DACB Saturation (V)' */
  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
    /* S-Function (hil_write_analog_block): '<S3>/HIL Write Analog' */

    /* S-Function Block: ExperimentRobustControl/SRV02+2DIP-E/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      result = hil_write_analog(ExperimentRobustControl_DW.HILInitialize_Card,
        ExperimentRobustControl_P.HILWriteAnalog_channels, 2,
        &ExperimentRobustControl_B.DACBSaturationV[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
      }
    }
  }

  /* Gain: '<S21>/Amplifier  Gain (V//V)' */
  y = -ExperimentRobustControl_P.K_AMP;
  ExperimentRobustControl_B.AmplifierGainVV[0] = y *
    ExperimentRobustControl_B.DACBSaturationV[0];
  ExperimentRobustControl_B.AmplifierGainVV[1] = y *
    ExperimentRobustControl_B.DACBSaturationV[1];

  /* Clock: '<S3>/Clock' */
  ExperimentRobustControl_B.Clock = ExperimentRobustControl_M->Timing.t[0];

  /* RelationalOperator: '<S22>/Compare' incorporates:
   *  Constant: '<S22>/Constant'
   */
  ExperimentRobustControl_B.Compare = (uint8_T)(ExperimentRobustControl_B.Clock <=
    ExperimentRobustControl_P.CompareToConstant_const);

  /* Outputs for Enabled SubSystem: '<S3>/Reset Encoders' incorporates:
   *  EnablePort: '<S25>/Enable'
   */
  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
    if (ExperimentRobustControl_B.Compare > 0) {
      if (!ExperimentRobustControl_DW.ResetEncoders_MODE) {
        ExperimentRobustControl_DW.ResetEncoders_MODE = true;
      }
    } else {
      if (ExperimentRobustControl_DW.ResetEncoders_MODE) {
        ExperimentRobustControl_DW.ResetEncoders_MODE = false;
      }
    }
  }

  if (ExperimentRobustControl_DW.ResetEncoders_MODE) {
    if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
      /* S-Function (hil_set_encoder_counts_block): '<S25>/HIL Set Encoder Counts' */

      /* S-Function Block: ExperimentRobustControl/SRV02+2DIP-E/Reset Encoders/HIL Set Encoder Counts (hil_set_encoder_counts_block) */
      {
        t_error result;
        ExperimentRobustControl_DW.HILSetEncoderCounts_Buffer[0] =
          ExperimentRobustControl_P.HILSetEncoderCounts_counts;
        ExperimentRobustControl_DW.HILSetEncoderCounts_Buffer[1] =
          ExperimentRobustControl_P.HILSetEncoderCounts_counts;
        ExperimentRobustControl_DW.HILSetEncoderCounts_Buffer[2] =
          ExperimentRobustControl_P.HILSetEncoderCounts_counts;
        ExperimentRobustControl_DW.HILSetEncoderCounts_Buffer[3] =
          ExperimentRobustControl_P.HILSetEncoderCounts_counts;
        result = hil_set_encoder_counts
          (ExperimentRobustControl_DW.HILInitialize_Card,
           ExperimentRobustControl_P.HILSetEncoderCounts_channels, 4,
           &ExperimentRobustControl_DW.HILSetEncoderCounts_Buffer[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        }
      }
    }

    if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
      srUpdateBC(ExperimentRobustControl_DW.ResetEncoders_SubsysRanBC);
    }
  }

  /* End of Outputs for SubSystem: '<S3>/Reset Encoders' */
  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
    /* S-Function (hil_write_digital_block): '<S3>/HIL Write Digital' */

    /* S-Function Block: ExperimentRobustControl/SRV02+2DIP-E/HIL Write Digital (hil_write_digital_block) */
    {
      t_error result;
      ExperimentRobustControl_DW.HILWriteDigital_Buffer[0] =
        (ExperimentRobustControl_P.EnableVoltPAQX2X4_Value[0] != 0);
      ExperimentRobustControl_DW.HILWriteDigital_Buffer[1] =
        (ExperimentRobustControl_P.EnableVoltPAQX2X4_Value[1] != 0);
      ExperimentRobustControl_DW.HILWriteDigital_Buffer[2] =
        (ExperimentRobustControl_P.EnableVoltPAQX2X4_Value[2] != 0);
      ExperimentRobustControl_DW.HILWriteDigital_Buffer[3] =
        (ExperimentRobustControl_P.EnableVoltPAQX2X4_Value[3] != 0);
      result = hil_write_digital(ExperimentRobustControl_DW.HILInitialize_Card,
        ExperimentRobustControl_P.HILWriteDigital_channels, 4,
        &ExperimentRobustControl_DW.HILWriteDigital_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
      }
    }

    /* Gain: '<Root>/Gain2' */
    ExperimentRobustControl_B.Gain2[0] = ExperimentRobustControl_P.Gain2_Gain *
      ExperimentRobustControl_B.EncoderCalibrationradcount;
    ExperimentRobustControl_B.Gain2[1] = ExperimentRobustControl_P.Gain2_Gain *
      ExperimentRobustControl_B.uDIPEncoderCalibrationradcount;

    /* Gain: '<Root>/Gain' */
    ExperimentRobustControl_B.Gain[0] = ExperimentRobustControl_P.Gain_Gain *
      ExperimentRobustControl_B.EncoderCalibrationYradcount;
    ExperimentRobustControl_B.Gain[1] = ExperimentRobustControl_P.Gain_Gain *
      ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount;
  }

  /* SignalConversion: '<S1>/TmpSignal ConversionAtState-SpaceInport1' */
  ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[0] =
    ExperimentRobustControl_B.EncoderCalibrationradcount;
  ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[1] =
    ExperimentRobustControl_B.uDIPEncoderCalibrationradcount;
  ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[2] =
    ExperimentRobustControl_B.EncoderCalibrationYradcount;
  ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[3] =
    ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount;

  /* Gain: '<S8>/Gain' */
  for (i = 0; i < 5; i++) {
    ExperimentRobustControl_B.Gain_l[i] = ExperimentRobustControl_P.Gain_Gain_d *
      ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[i];
  }

  /* End of Gain: '<S8>/Gain' */
  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
  }

  /* Gain: '<S4>/-sysc.c_x' */
  for (i = 0; i < 2; i++) {
    ExperimentRobustControl_B.syscc_x[i] = 0.0;
    ExperimentRobustControl_B.syscc_x[i] +=
      ExperimentRobustControl_P.syscc_x_Gain[i] *
      ExperimentRobustControl_B.Integrator[0];
    ExperimentRobustControl_B.syscc_x[i] +=
      ExperimentRobustControl_P.syscc_x_Gain[i + 2] *
      ExperimentRobustControl_B.Integrator[1];
    ExperimentRobustControl_B.syscc_x[i] +=
      ExperimentRobustControl_P.syscc_x_Gain[i + 4] *
      ExperimentRobustControl_B.Integrator[2];
    ExperimentRobustControl_B.syscc_x[i] +=
      ExperimentRobustControl_P.syscc_x_Gain[i + 6] *
      ExperimentRobustControl_B.Integrator[3];
  }

  /* End of Gain: '<S4>/-sysc.c_x' */

  /* Sum: '<S4>/Sum1' */
  ExperimentRobustControl_B.Sum1[0] =
    ExperimentRobustControl_B.EncoderCalibrationradcount -
    ExperimentRobustControl_B.syscc_x[0];
  ExperimentRobustControl_B.Sum1[1] =
    ExperimentRobustControl_B.uDIPEncoderCalibrationradcount -
    ExperimentRobustControl_B.syscc_x[1];
  for (i = 0; i < 4; i++) {
    /* Gain: '<S4>/Kalman Filter Gain' */
    ExperimentRobustControl_B.KalmanFilterGain[i] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain[i] +=
      ExperimentRobustControl_P.L_CT_labb[i] * ExperimentRobustControl_B.Sum1[0];
    ExperimentRobustControl_B.KalmanFilterGain[i] +=
      ExperimentRobustControl_P.L_CT_labb[i + 4] *
      ExperimentRobustControl_B.Sum1[1];

    /* Gain: '<S4>/sysc.b_x' */
    ExperimentRobustControl_B.syscb_x[i] =
      ExperimentRobustControl_P.syscb_x_Gain[i] *
      ExperimentRobustControl_B.ManualSwitch2[0];

    /* Gain: '<S4>/sysc.a_x' */
    ExperimentRobustControl_B.sysca_x[i] = 0.0;
    ExperimentRobustControl_B.sysca_x[i] +=
      ExperimentRobustControl_P.sysca_x_Gain[i] *
      ExperimentRobustControl_B.Integrator[0];
    ExperimentRobustControl_B.sysca_x[i] +=
      ExperimentRobustControl_P.sysca_x_Gain[i + 4] *
      ExperimentRobustControl_B.Integrator[1];
    ExperimentRobustControl_B.sysca_x[i] +=
      ExperimentRobustControl_P.sysca_x_Gain[i + 8] *
      ExperimentRobustControl_B.Integrator[2];
    ExperimentRobustControl_B.sysca_x[i] +=
      ExperimentRobustControl_P.sysca_x_Gain[i + 12] *
      ExperimentRobustControl_B.Integrator[3];

    /* Sum: '<S4>/Sum' */
    ExperimentRobustControl_B.Sum[i] =
      (ExperimentRobustControl_B.KalmanFilterGain[i] +
       ExperimentRobustControl_B.syscb_x[i]) +
      ExperimentRobustControl_B.sysca_x[i];
  }

  /* Gain: '<S5>/-sysd.c_x' */
  for (i = 0; i < 2; i++) {
    ExperimentRobustControl_B.sysdc_x[i] = 0.0;
    ExperimentRobustControl_B.sysdc_x[i] +=
      ExperimentRobustControl_P.sysdc_x_Gain[i] *
      ExperimentRobustControl_B.Integrator_k[0];
    ExperimentRobustControl_B.sysdc_x[i] +=
      ExperimentRobustControl_P.sysdc_x_Gain[i + 2] *
      ExperimentRobustControl_B.Integrator_k[1];
    ExperimentRobustControl_B.sysdc_x[i] +=
      ExperimentRobustControl_P.sysdc_x_Gain[i + 4] *
      ExperimentRobustControl_B.Integrator_k[2];
    ExperimentRobustControl_B.sysdc_x[i] +=
      ExperimentRobustControl_P.sysdc_x_Gain[i + 6] *
      ExperimentRobustControl_B.Integrator_k[3];
  }

  /* End of Gain: '<S5>/-sysd.c_x' */

  /* Sum: '<S5>/Sum1' */
  ExperimentRobustControl_B.Sum1_k[0] =
    ExperimentRobustControl_B.EncoderCalibrationYradcount -
    ExperimentRobustControl_B.sysdc_x[0];
  ExperimentRobustControl_B.Sum1_k[1] =
    ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount -
    ExperimentRobustControl_B.sysdc_x[1];
  for (i = 0; i < 4; i++) {
    /* Gain: '<S5>/Kalman Filter Gain' */
    ExperimentRobustControl_B.KalmanFilterGain_e[i] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain_e[i] +=
      ExperimentRobustControl_P.L_CT_labb[i] * ExperimentRobustControl_B.Sum1_k
      [0];
    ExperimentRobustControl_B.KalmanFilterGain_e[i] +=
      ExperimentRobustControl_P.L_CT_labb[i + 4] *
      ExperimentRobustControl_B.Sum1_k[1];

    /* Gain: '<S5>/sysd.b_x' */
    ExperimentRobustControl_B.sysdb_x[i] =
      ExperimentRobustControl_P.sysdb_x_Gain[i] *
      ExperimentRobustControl_B.ManualSwitch2[1];

    /* Gain: '<S5>/sysd.a_x' */
    ExperimentRobustControl_B.sysda_x[i] = 0.0;
    ExperimentRobustControl_B.sysda_x[i] +=
      ExperimentRobustControl_P.sysda_x_Gain[i] *
      ExperimentRobustControl_B.Integrator_k[0];
    ExperimentRobustControl_B.sysda_x[i] +=
      ExperimentRobustControl_P.sysda_x_Gain[i + 4] *
      ExperimentRobustControl_B.Integrator_k[1];
    ExperimentRobustControl_B.sysda_x[i] +=
      ExperimentRobustControl_P.sysda_x_Gain[i + 8] *
      ExperimentRobustControl_B.Integrator_k[2];
    ExperimentRobustControl_B.sysda_x[i] +=
      ExperimentRobustControl_P.sysda_x_Gain[i + 12] *
      ExperimentRobustControl_B.Integrator_k[3];

    /* Sum: '<S5>/Sum' */
    ExperimentRobustControl_B.Sum_j[i] =
      (ExperimentRobustControl_B.KalmanFilterGain_e[i] +
       ExperimentRobustControl_B.sysdb_x[i]) +
      ExperimentRobustControl_B.sysda_x[i];
  }

  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
    /* Outputs for Atomic SubSystem: '<S6>/Pendulum Angle Watchdog: X' */
    Experime_PendulumAngleWatchdogX(ExperimentRobustControl_M,
      ExperimentRobustControl_B.uDIPEncoderCalibrationradcount,
      &ExperimentRobustControl_B.PendulumAngleWatchdogX_d,
      (P_PendulumAngleWatchdogX_Expe_T *)
      &ExperimentRobustControl_P.PendulumAngleWatchdogX_d);

    /* End of Outputs for SubSystem: '<S6>/Pendulum Angle Watchdog: X' */

    /* Outputs for Atomic SubSystem: '<S6>/SRV02 X Angle Watchdog' */
    ExperimentR_SRV02XAngleWatchdog(ExperimentRobustControl_M,
      ExperimentRobustControl_B.EncoderCalibrationradcount,
      &ExperimentRobustControl_B.SRV02XAngleWatchdog_n,
      (P_SRV02XAngleWatchdog_Experim_T *)
      &ExperimentRobustControl_P.SRV02XAngleWatchdog_n);

    /* End of Outputs for SubSystem: '<S6>/SRV02 X Angle Watchdog' */

    /* Outputs for Atomic SubSystem: '<S7>/Pendulum Angle Watchdog: Y' */
    Experime_PendulumAngleWatchdogY(ExperimentRobustControl_M,
      ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount,
      &ExperimentRobustControl_B.PendulumAngleWatchdogY_a,
      (P_PendulumAngleWatchdogY_Expe_T *)
      &ExperimentRobustControl_P.PendulumAngleWatchdogY_a);

    /* End of Outputs for SubSystem: '<S7>/Pendulum Angle Watchdog: Y' */

    /* Outputs for Atomic SubSystem: '<S7>/SRV02 Y Angle Watchdog' */
    ExperimentR_SRV02YAngleWatchdog(ExperimentRobustControl_M,
      ExperimentRobustControl_B.EncoderCalibrationYradcount,
      &ExperimentRobustControl_B.SRV02YAngleWatchdog_l,
      (P_SRV02YAngleWatchdog_Experim_T *)
      &ExperimentRobustControl_P.SRV02YAngleWatchdog_l);

    /* End of Outputs for SubSystem: '<S7>/SRV02 Y Angle Watchdog' */
  }
}

/* Model update function */
void ExperimentRobustControl_update(void)
{
  if (rtmIsMajorTimeStep(ExperimentRobustControl_M)) {
    rt_ertODEUpdateContinuousStates(&ExperimentRobustControl_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++ExperimentRobustControl_M->Timing.clockTick0)) {
    ++ExperimentRobustControl_M->Timing.clockTickH0;
  }

  ExperimentRobustControl_M->Timing.t[0] = rtsiGetSolverStopTime
    (&ExperimentRobustControl_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++ExperimentRobustControl_M->Timing.clockTick1)) {
      ++ExperimentRobustControl_M->Timing.clockTickH1;
    }

    ExperimentRobustControl_M->Timing.t[1] =
      ExperimentRobustControl_M->Timing.clockTick1 *
      ExperimentRobustControl_M->Timing.stepSize1 +
      ExperimentRobustControl_M->Timing.clockTickH1 *
      ExperimentRobustControl_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void ExperimentRobustControl_derivatives(void)
{
  int_T is;
  int_T ci;
  XDot_ExperimentRobustControl_T *_rtXdot;
  _rtXdot = ((XDot_ExperimentRobustControl_T *)
             ExperimentRobustControl_M->derivs);

  /* Derivatives for StateSpace: '<S1>/State-Space' */
  for (is = 0; is < 54; is++) {
    _rtXdot->StateSpace_CSTATE[is] = 0.0;
    for (ci = 0; ci < 54; ci++) {
      _rtXdot->StateSpace_CSTATE[is] += ExperimentRobustControl_P.ah[ci * 54 +
        is] * ExperimentRobustControl_X.StateSpace_CSTATE[ci];
    }

    _rtXdot->StateSpace_CSTATE[is] += ExperimentRobustControl_P.bh[is] *
      ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[0];
    _rtXdot->StateSpace_CSTATE[is] += ExperimentRobustControl_P.bh[54 + is] *
      ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[1];
    _rtXdot->StateSpace_CSTATE[is] += ExperimentRobustControl_P.bh[108 + is] *
      ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[2];
    _rtXdot->StateSpace_CSTATE[is] += ExperimentRobustControl_P.bh[162 + is] *
      ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[3];
  }

  /* End of Derivatives for StateSpace: '<S1>/State-Space' */

  /* Derivatives for Integrator: '<S6>/Integrator' */
  _rtXdot->Integrator_CSTATE_c =
    ExperimentRobustControl_B.EncoderCalibrationradcount;

  /* Derivatives for Integrator: '<S4>/Integrator' */
  _rtXdot->Integrator_CSTATE[0] = ExperimentRobustControl_B.Sum[0];

  /* Derivatives for Integrator: '<S5>/Integrator' */
  _rtXdot->Integrator_CSTATE_n[0] = ExperimentRobustControl_B.Sum_j[0];

  /* Derivatives for Integrator: '<S4>/Integrator' */
  _rtXdot->Integrator_CSTATE[1] = ExperimentRobustControl_B.Sum[1];

  /* Derivatives for Integrator: '<S5>/Integrator' */
  _rtXdot->Integrator_CSTATE_n[1] = ExperimentRobustControl_B.Sum_j[1];

  /* Derivatives for Integrator: '<S4>/Integrator' */
  _rtXdot->Integrator_CSTATE[2] = ExperimentRobustControl_B.Sum[2];

  /* Derivatives for Integrator: '<S5>/Integrator' */
  _rtXdot->Integrator_CSTATE_n[2] = ExperimentRobustControl_B.Sum_j[2];

  /* Derivatives for Integrator: '<S4>/Integrator' */
  _rtXdot->Integrator_CSTATE[3] = ExperimentRobustControl_B.Sum[3];

  /* Derivatives for Integrator: '<S5>/Integrator' */
  _rtXdot->Integrator_CSTATE_n[3] = ExperimentRobustControl_B.Sum_j[3];

  /* Derivatives for Integrator: '<S7>/Integrator' */
  _rtXdot->Integrator_CSTATE_a =
    ExperimentRobustControl_B.EncoderCalibrationYradcount;
}

/* Model initialize function */
void ExperimentRobustControl_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<S3>/HIL Initialize' */

  /* S-Function Block: ExperimentRobustControl/SRV02+2DIP-E/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0",
                      &ExperimentRobustControl_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (ExperimentRobustControl_DW.HILInitialize_Card, "update_rate=normal", 19);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(ExperimentRobustControl_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
      return;
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_analog_input_ &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_analog_inpu_o &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &ExperimentRobustControl_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            ExperimentRobustControl_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &ExperimentRobustControl_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            ExperimentRobustControl_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_analog_input_chan, 8U,
         &ExperimentRobustControl_DW.HILInitialize_AIMinimums[0],
         &ExperimentRobustControl_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_analog_output &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_analog_outp_d &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums =
          &ExperimentRobustControl_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            ExperimentRobustControl_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums =
          &ExperimentRobustControl_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            ExperimentRobustControl_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_analog_output_cha, 8U,
         &ExperimentRobustControl_DW.HILInitialize_AOMinimums[0],
         &ExperimentRobustControl_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_analog_outp_j &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_analog_outp_m &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &ExperimentRobustControl_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            ExperimentRobustControl_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(ExperimentRobustControl_DW.HILInitialize_Card,
        ExperimentRobustControl_P.HILInitialize_analog_output_cha, 8U,
        &ExperimentRobustControl_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if (ExperimentRobustControl_P.HILInitialize_set_analog_outp_k) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &ExperimentRobustControl_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            ExperimentRobustControl_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_analog_output_cha, 8U,
         &ExperimentRobustControl_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    result = hil_set_digital_directions
      (ExperimentRobustControl_DW.HILInitialize_Card, NULL, 0U,
       ExperimentRobustControl_P.HILInitialize_digital_output_ch, 8U);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
      return;
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_digital_out_a &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_digital_out_f &&
         is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits = &ExperimentRobustControl_DW.HILInitialize_DOBits
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] =
            ExperimentRobustControl_P.HILInitialize_initial_digital_o;
        }
      }

      result = hil_write_digital(ExperimentRobustControl_DW.HILInitialize_Card,
        ExperimentRobustControl_P.HILInitialize_digital_output_ch, 8U,
        (t_boolean *) &ExperimentRobustControl_DW.HILInitialize_DOBits[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if (ExperimentRobustControl_P.HILInitialize_set_digital_ou_j0) {
      {
        int_T i1;
        int32_T *dw_DOStates =
          &ExperimentRobustControl_DW.HILInitialize_DOStates[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOStates[i1] =
            ExperimentRobustControl_P.HILInitialize_watchdog_digital_;
        }
      }

      result = hil_watchdog_set_digital_expiration_state
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_digital_output_ch, 8U, (const
          t_digital_state *) &ExperimentRobustControl_DW.HILInitialize_DOStates
         [0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_encoder_param &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_encoder_par_m &&
         is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &ExperimentRobustControl_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] =
            ExperimentRobustControl_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &ExperimentRobustControl_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_encoder_count &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_encoder_cou_e &&
         is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &ExperimentRobustControl_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            ExperimentRobustControl_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_encoder_channels, 8U,
         &ExperimentRobustControl_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_pwm_params_at &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_pwm_params__m &&
         is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &ExperimentRobustControl_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            ExperimentRobustControl_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(ExperimentRobustControl_DW.HILInitialize_Card,
        ExperimentRobustControl_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &ExperimentRobustControl_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          ExperimentRobustControl_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &ExperimentRobustControl_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            ExperimentRobustControl_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            ExperimentRobustControl_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = ExperimentRobustControl_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            ExperimentRobustControl_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            ExperimentRobustControl_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              ExperimentRobustControl_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency
          (ExperimentRobustControl_DW.HILInitialize_Card,
           &ExperimentRobustControl_DW.HILInitialize_POSortedChans[0],
           num_duty_cycle_modes,
           &ExperimentRobustControl_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle
          (ExperimentRobustControl_DW.HILInitialize_Card,
           &ExperimentRobustControl_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
           num_frequency_modes,
           &ExperimentRobustControl_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &ExperimentRobustControl_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            ExperimentRobustControl_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &ExperimentRobustControl_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] =
            ExperimentRobustControl_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &ExperimentRobustControl_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] =
            ExperimentRobustControl_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_pwm_channels, 8U,
         (t_pwm_configuration *)
         &ExperimentRobustControl_DW.HILInitialize_POModeValues[0],
         (t_pwm_alignment *)
         &ExperimentRobustControl_DW.HILInitialize_POAlignValues[0],
         (t_pwm_polarity *)
         &ExperimentRobustControl_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &ExperimentRobustControl_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            ExperimentRobustControl_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues =
          &ExperimentRobustControl_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            ExperimentRobustControl_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_pwm_channels, 8U,
         &ExperimentRobustControl_DW.HILInitialize_POSortedFreqs[0],
         &ExperimentRobustControl_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_a &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_l &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues =
          &ExperimentRobustControl_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            ExperimentRobustControl_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(ExperimentRobustControl_DW.HILInitialize_Card,
        ExperimentRobustControl_P.HILInitialize_pwm_channels, 8U,
        &ExperimentRobustControl_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }

    if (ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues =
          &ExperimentRobustControl_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            ExperimentRobustControl_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (ExperimentRobustControl_DW.HILInitialize_Card,
         ExperimentRobustControl_P.HILInitialize_pwm_channels, 8U,
         &ExperimentRobustControl_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S3>/HIL Read Encoder Timebase' */

  /* S-Function Block: ExperimentRobustControl/SRV02+2DIP-E/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (ExperimentRobustControl_DW.HILInitialize_Card,
       ExperimentRobustControl_P.HILReadEncoderTimebase_samples_,
       ExperimentRobustControl_P.HILReadEncoderTimebase_channels, 4,
       &ExperimentRobustControl_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
    }
  }

  /* Start for Enabled SubSystem: '<S3>/Reset Encoders' */
  ExperimentRobustControl_DW.ResetEncoders_MODE = false;

  /* End of Start for SubSystem: '<S3>/Reset Encoders' */
  {
    int_T is;

    /* InitializeConditions for StateSpace: '<S1>/State-Space' */
    for (is = 0; is < 54; is++) {
      ExperimentRobustControl_X.StateSpace_CSTATE[is] =
        ExperimentRobustControl_P.StateSpace_X0;
    }

    /* End of InitializeConditions for StateSpace: '<S1>/State-Space' */

    /* InitializeConditions for Integrator: '<S6>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE_c =
      ExperimentRobustControl_P.Integrator_IC_b;

    /* InitializeConditions for Integrator: '<S4>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE[0] =
      ExperimentRobustControl_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE_n[0] =
      ExperimentRobustControl_P.Integrator_IC_p;

    /* InitializeConditions for Integrator: '<S4>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE[1] =
      ExperimentRobustControl_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE_n[1] =
      ExperimentRobustControl_P.Integrator_IC_p;

    /* InitializeConditions for Integrator: '<S4>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE[2] =
      ExperimentRobustControl_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE_n[2] =
      ExperimentRobustControl_P.Integrator_IC_p;

    /* InitializeConditions for Integrator: '<S4>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE[3] =
      ExperimentRobustControl_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S5>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE_n[3] =
      ExperimentRobustControl_P.Integrator_IC_p;

    /* InitializeConditions for Integrator: '<S7>/Integrator' */
    ExperimentRobustControl_X.Integrator_CSTATE_a =
      ExperimentRobustControl_P.Integrator_IC_px;
  }
}

/* Model terminate function */
void ExperimentRobustControl_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S3>/HIL Initialize' */

  /* S-Function Block: ExperimentRobustControl/SRV02+2DIP-E/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_digital_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(ExperimentRobustControl_DW.HILInitialize_Card);
    hil_monitor_stop_all(ExperimentRobustControl_DW.HILInitialize_Card);
    is_switching = false;
    if ((ExperimentRobustControl_P.HILInitialize_set_analog_outp_n &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_analog_out_mi &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &ExperimentRobustControl_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            ExperimentRobustControl_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_digital_ou_bw &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_digital_out_b &&
         is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits = &ExperimentRobustControl_DW.HILInitialize_DOBits
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] =
            ExperimentRobustControl_P.HILInitialize_final_digital_out;
        }
      }

      num_final_digital_outputs = 8U;
    }

    if ((ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_c &&
         !is_switching) ||
        (ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_n &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues =
          &ExperimentRobustControl_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            ExperimentRobustControl_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        || num_final_digital_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(ExperimentRobustControl_DW.HILInitialize_Card
                         ,
                         ExperimentRobustControl_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , ExperimentRobustControl_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         ,
                         ExperimentRobustControl_P.HILInitialize_digital_output_ch,
                         num_final_digital_outputs
                         , NULL, 0
                         , &ExperimentRobustControl_DW.HILInitialize_AOVoltages
                         [0]
                         , &ExperimentRobustControl_DW.HILInitialize_POValues[0]
                         , (t_boolean *)
                         &ExperimentRobustControl_DW.HILInitialize_DOBits[0]
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (ExperimentRobustControl_DW.HILInitialize_Card,
             ExperimentRobustControl_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &ExperimentRobustControl_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm
            (ExperimentRobustControl_DW.HILInitialize_Card,
             ExperimentRobustControl_P.HILInitialize_pwm_channels,
             num_final_pwm_outputs,
             &ExperimentRobustControl_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_digital_outputs > 0) {
          local_result = hil_write_digital
            (ExperimentRobustControl_DW.HILInitialize_Card,
             ExperimentRobustControl_P.HILInitialize_digital_output_ch,
             num_final_digital_outputs, (t_boolean *)
             &ExperimentRobustControl_DW.HILInitialize_DOBits[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(ExperimentRobustControl_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(ExperimentRobustControl_DW.HILInitialize_Card);
    hil_monitor_delete_all(ExperimentRobustControl_DW.HILInitialize_Card);
    hil_close(ExperimentRobustControl_DW.HILInitialize_Card);
    ExperimentRobustControl_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  ExperimentRobustControl_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  ExperimentRobustControl_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  ExperimentRobustControl_initialize();
}

void MdlTerminate(void)
{
  ExperimentRobustControl_terminate();
}

/* Registration function */
RT_MODEL_ExperimentRobustCont_T *ExperimentRobustControl(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)ExperimentRobustControl_M, 0,
                sizeof(RT_MODEL_ExperimentRobustCont_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&ExperimentRobustControl_M->solverInfo,
                          &ExperimentRobustControl_M->Timing.simTimeStep);
    rtsiSetTPtr(&ExperimentRobustControl_M->solverInfo, &rtmGetTPtr
                (ExperimentRobustControl_M));
    rtsiSetStepSizePtr(&ExperimentRobustControl_M->solverInfo,
                       &ExperimentRobustControl_M->Timing.stepSize0);
    rtsiSetdXPtr(&ExperimentRobustControl_M->solverInfo,
                 &ExperimentRobustControl_M->derivs);
    rtsiSetContStatesPtr(&ExperimentRobustControl_M->solverInfo, (real_T **)
                         &ExperimentRobustControl_M->contStates);
    rtsiSetNumContStatesPtr(&ExperimentRobustControl_M->solverInfo,
      &ExperimentRobustControl_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&ExperimentRobustControl_M->solverInfo,
      &ExperimentRobustControl_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&ExperimentRobustControl_M->solverInfo,
      &ExperimentRobustControl_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&ExperimentRobustControl_M->solverInfo,
      &ExperimentRobustControl_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&ExperimentRobustControl_M->solverInfo,
                          (&rtmGetErrorStatus(ExperimentRobustControl_M)));
    rtsiSetRTModelPtr(&ExperimentRobustControl_M->solverInfo,
                      ExperimentRobustControl_M);
  }

  rtsiSetSimTimeStep(&ExperimentRobustControl_M->solverInfo, MAJOR_TIME_STEP);
  ExperimentRobustControl_M->intgData.f[0] = ExperimentRobustControl_M->odeF[0];
  ExperimentRobustControl_M->contStates = ((real_T *) &ExperimentRobustControl_X);
  rtsiSetSolverData(&ExperimentRobustControl_M->solverInfo, (void *)
                    &ExperimentRobustControl_M->intgData);
  rtsiSetSolverName(&ExperimentRobustControl_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = ExperimentRobustControl_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    ExperimentRobustControl_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    ExperimentRobustControl_M->Timing.sampleTimes =
      (&ExperimentRobustControl_M->Timing.sampleTimesArray[0]);
    ExperimentRobustControl_M->Timing.offsetTimes =
      (&ExperimentRobustControl_M->Timing.offsetTimesArray[0]);

    /* task periods */
    ExperimentRobustControl_M->Timing.sampleTimes[0] = (0.0);
    ExperimentRobustControl_M->Timing.sampleTimes[1] = (0.001);

    /* task offsets */
    ExperimentRobustControl_M->Timing.offsetTimes[0] = (0.0);
    ExperimentRobustControl_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(ExperimentRobustControl_M,
             &ExperimentRobustControl_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = ExperimentRobustControl_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    ExperimentRobustControl_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(ExperimentRobustControl_M, -1);
  ExperimentRobustControl_M->Timing.stepSize0 = 0.001;
  ExperimentRobustControl_M->Timing.stepSize1 = 0.001;

  /* External mode info */
  ExperimentRobustControl_M->Sizes.checksums[0] = (3116169259U);
  ExperimentRobustControl_M->Sizes.checksums[1] = (2044929189U);
  ExperimentRobustControl_M->Sizes.checksums[2] = (3951391620U);
  ExperimentRobustControl_M->Sizes.checksums[3] = (600318527U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[11];
    ExperimentRobustControl_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    systemRan[3] = &rtAlwaysEnabled;
    systemRan[4] = &rtAlwaysEnabled;
    systemRan[5] = &rtAlwaysEnabled;
    systemRan[6] = &rtAlwaysEnabled;
    systemRan[7] = &rtAlwaysEnabled;
    systemRan[8] = &rtAlwaysEnabled;
    systemRan[9] = &rtAlwaysEnabled;
    systemRan[10] = (sysRanDType *)
      &ExperimentRobustControl_DW.ResetEncoders_SubsysRanBC;
    rteiSetModelMappingInfoPtr(ExperimentRobustControl_M->extModeInfo,
      &ExperimentRobustControl_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(ExperimentRobustControl_M->extModeInfo,
                        ExperimentRobustControl_M->Sizes.checksums);
    rteiSetTPtr(ExperimentRobustControl_M->extModeInfo, rtmGetTPtr
                (ExperimentRobustControl_M));
  }

  ExperimentRobustControl_M->solverInfoPtr =
    (&ExperimentRobustControl_M->solverInfo);
  ExperimentRobustControl_M->Timing.stepSize = (0.001);
  rtsiSetFixedStepSize(&ExperimentRobustControl_M->solverInfo, 0.001);
  rtsiSetSolverMode(&ExperimentRobustControl_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  ExperimentRobustControl_M->blockIO = ((void *) &ExperimentRobustControl_B);
  (void) memset(((void *) &ExperimentRobustControl_B), 0,
                sizeof(B_ExperimentRobustControl_T));

  {
    int32_T i;
    for (i = 0; i < 5; i++) {
      ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn[i] = 0.0;
    }

    for (i = 0; i < 5; i++) {
      ExperimentRobustControl_B.Gain_l[i] = 0.0;
    }

    for (i = 0; i < 5; i++) {
      ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn[i] = 0.0;
    }

    ExperimentRobustControl_B.HILReadEncoderTimebase_o1 = 0.0;
    ExperimentRobustControl_B.HILReadEncoderTimebase_o2 = 0.0;
    ExperimentRobustControl_B.HILReadEncoderTimebase_o3 = 0.0;
    ExperimentRobustControl_B.HILReadEncoderTimebase_o4 = 0.0;
    ExperimentRobustControl_B.uDIPEncoderCalibrationradcount = 0.0;
    ExperimentRobustControl_B.EncoderCalibrationradcount = 0.0;
    ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount = 0.0;
    ExperimentRobustControl_B.EncoderCalibrationYradcount = 0.0;
    ExperimentRobustControl_B.StateSpace[0] = 0.0;
    ExperimentRobustControl_B.StateSpace[1] = 0.0;
    ExperimentRobustControl_B.Integrator[0] = 0.0;
    ExperimentRobustControl_B.Integrator[1] = 0.0;
    ExperimentRobustControl_B.Integrator[2] = 0.0;
    ExperimentRobustControl_B.Integrator[3] = 0.0;
    ExperimentRobustControl_B.Integrator_m = 0.0;
    ExperimentRobustControl_B.Integrator_k[0] = 0.0;
    ExperimentRobustControl_B.Integrator_k[1] = 0.0;
    ExperimentRobustControl_B.Integrator_k[2] = 0.0;
    ExperimentRobustControl_B.Integrator_k[3] = 0.0;
    ExperimentRobustControl_B.Integrator_c = 0.0;
    ExperimentRobustControl_B.ManualSwitch2[0] = 0.0;
    ExperimentRobustControl_B.ManualSwitch2[1] = 0.0;
    ExperimentRobustControl_B.Gain1[0] = 0.0;
    ExperimentRobustControl_B.Gain1[1] = 0.0;
    ExperimentRobustControl_B.DirectionConventionRightHandsys[0] = 0.0;
    ExperimentRobustControl_B.DirectionConventionRightHandsys[1] = 0.0;
    ExperimentRobustControl_B.AmplifierSaturationV[0] = 0.0;
    ExperimentRobustControl_B.AmplifierSaturationV[1] = 0.0;
    ExperimentRobustControl_B.InverseAmplifierGainVV[0] = 0.0;
    ExperimentRobustControl_B.InverseAmplifierGainVV[1] = 0.0;
    ExperimentRobustControl_B.DACBSaturationV[0] = 0.0;
    ExperimentRobustControl_B.DACBSaturationV[1] = 0.0;
    ExperimentRobustControl_B.AmplifierGainVV[0] = 0.0;
    ExperimentRobustControl_B.AmplifierGainVV[1] = 0.0;
    ExperimentRobustControl_B.Clock = 0.0;
    ExperimentRobustControl_B.Gain2[0] = 0.0;
    ExperimentRobustControl_B.Gain2[1] = 0.0;
    ExperimentRobustControl_B.Gain[0] = 0.0;
    ExperimentRobustControl_B.Gain[1] = 0.0;
    ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[0] = 0.0;
    ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[1] = 0.0;
    ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[2] = 0.0;
    ExperimentRobustControl_B.TmpSignalConversionAtStateSpace[3] = 0.0;
    ExperimentRobustControl_B.syscc_x[0] = 0.0;
    ExperimentRobustControl_B.syscc_x[1] = 0.0;
    ExperimentRobustControl_B.Sum1[0] = 0.0;
    ExperimentRobustControl_B.Sum1[1] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain[0] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain[1] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain[2] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain[3] = 0.0;
    ExperimentRobustControl_B.syscb_x[0] = 0.0;
    ExperimentRobustControl_B.syscb_x[1] = 0.0;
    ExperimentRobustControl_B.syscb_x[2] = 0.0;
    ExperimentRobustControl_B.syscb_x[3] = 0.0;
    ExperimentRobustControl_B.sysca_x[0] = 0.0;
    ExperimentRobustControl_B.sysca_x[1] = 0.0;
    ExperimentRobustControl_B.sysca_x[2] = 0.0;
    ExperimentRobustControl_B.sysca_x[3] = 0.0;
    ExperimentRobustControl_B.Sum[0] = 0.0;
    ExperimentRobustControl_B.Sum[1] = 0.0;
    ExperimentRobustControl_B.Sum[2] = 0.0;
    ExperimentRobustControl_B.Sum[3] = 0.0;
    ExperimentRobustControl_B.sysdc_x[0] = 0.0;
    ExperimentRobustControl_B.sysdc_x[1] = 0.0;
    ExperimentRobustControl_B.Sum1_k[0] = 0.0;
    ExperimentRobustControl_B.Sum1_k[1] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain_e[0] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain_e[1] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain_e[2] = 0.0;
    ExperimentRobustControl_B.KalmanFilterGain_e[3] = 0.0;
    ExperimentRobustControl_B.sysdb_x[0] = 0.0;
    ExperimentRobustControl_B.sysdb_x[1] = 0.0;
    ExperimentRobustControl_B.sysdb_x[2] = 0.0;
    ExperimentRobustControl_B.sysdb_x[3] = 0.0;
    ExperimentRobustControl_B.sysda_x[0] = 0.0;
    ExperimentRobustControl_B.sysda_x[1] = 0.0;
    ExperimentRobustControl_B.sysda_x[2] = 0.0;
    ExperimentRobustControl_B.sysda_x[3] = 0.0;
    ExperimentRobustControl_B.Sum_j[0] = 0.0;
    ExperimentRobustControl_B.Sum_j[1] = 0.0;
    ExperimentRobustControl_B.Sum_j[2] = 0.0;
    ExperimentRobustControl_B.Sum_j[3] = 0.0;
    ExperimentRobustControl_B.LQRGainX = 0.0;
    ExperimentRobustControl_B.LQRGainY = 0.0;
    ExperimentRobustControl_B.SRV02YAngleWatchdog.Abs = 0.0;
    ExperimentRobustControl_B.SRV02YAngleWatchdog.CheckMaximumAllowedAngle = 0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogY.Abs = 0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogY.CheckMaximumAllowedAngle =
      0.0;
    ExperimentRobustControl_B.SRV02XAngleWatchdog.Abs = 0.0;
    ExperimentRobustControl_B.SRV02XAngleWatchdog.CheckMaximumAllowedAngle = 0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogX.Abs = 0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogX.CheckMaximumAllowedAngle =
      0.0;
    ExperimentRobustControl_B.SRV02YAngleWatchdog_l.Abs = 0.0;
    ExperimentRobustControl_B.SRV02YAngleWatchdog_l.CheckMaximumAllowedAngle =
      0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogY_a.Abs = 0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogY_a.CheckMaximumAllowedAngle =
      0.0;
    ExperimentRobustControl_B.SRV02XAngleWatchdog_n.Abs = 0.0;
    ExperimentRobustControl_B.SRV02XAngleWatchdog_n.CheckMaximumAllowedAngle =
      0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogX_d.Abs = 0.0;
    ExperimentRobustControl_B.PendulumAngleWatchdogX_d.CheckMaximumAllowedAngle =
      0.0;
  }

  /* parameters */
  ExperimentRobustControl_M->defaultParam = ((real_T *)
    &ExperimentRobustControl_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &ExperimentRobustControl_X;
    ExperimentRobustControl_M->contStates = (x);
    (void) memset((void *)&ExperimentRobustControl_X, 0,
                  sizeof(X_ExperimentRobustControl_T));
  }

  /* states (dwork) */
  ExperimentRobustControl_M->dwork = ((void *) &ExperimentRobustControl_DW);
  (void) memset((void *)&ExperimentRobustControl_DW, 0,
                sizeof(DW_ExperimentRobustControl_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      ExperimentRobustControl_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    ExperimentRobustControl_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 17;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  ExperimentRobustControl_M->Sizes.numContStates = (64);/* Number of continuous states */
  ExperimentRobustControl_M->Sizes.numPeriodicContStates = (0);/* Number of periodic continuous states */
  ExperimentRobustControl_M->Sizes.numY = (0);/* Number of model outputs */
  ExperimentRobustControl_M->Sizes.numU = (0);/* Number of model inputs */
  ExperimentRobustControl_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  ExperimentRobustControl_M->Sizes.numSampTimes = (2);/* Number of sample times */
  ExperimentRobustControl_M->Sizes.numBlocks = (85);/* Number of blocks */
  ExperimentRobustControl_M->Sizes.numBlockIO = (66);/* Number of block outputs */
  ExperimentRobustControl_M->Sizes.numBlockPrms = (3472);/* Sum of parameter "widths" */
  return ExperimentRobustControl_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
