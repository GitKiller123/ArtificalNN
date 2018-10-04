/*
 * ExperimentRobustControl_dt.h
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

#include "ext_types.h"

/* data type size table */
static uint_T rtDataTypeSizes[] = {
  sizeof(real_T),
  sizeof(real32_T),
  sizeof(int8_T),
  sizeof(uint8_T),
  sizeof(int16_T),
  sizeof(uint16_T),
  sizeof(int32_T),
  sizeof(uint32_T),
  sizeof(boolean_T),
  sizeof(fcn_call_T),
  sizeof(int_T),
  sizeof(pointer_T),
  sizeof(action_T),
  2*sizeof(uint32_T),
  sizeof(t_card),
  sizeof(t_task),
  sizeof(t_boolean)
};

/* data type name table */
static const char_T * rtDataTypeNames[] = {
  "real_T",
  "real32_T",
  "int8_T",
  "uint8_T",
  "int16_T",
  "uint16_T",
  "int32_T",
  "uint32_T",
  "boolean_T",
  "fcn_call_T",
  "int_T",
  "pointer_T",
  "action_T",
  "timer_uint32_pair_T",
  "t_card",
  "t_task",
  "t_boolean"
};

/* data type transitions for block I/O structure */
static DataTypeTransition rtBTransitions[] = {
  { (char_T *)(&ExperimentRobustControl_B.HILReadEncoderTimebase_o1), 0, 0, 100
  },

  { (char_T *)(&ExperimentRobustControl_B.Compare), 3, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02YAngleWatchdog.Abs), 0, 0, 2 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02YAngleWatchdog.Compare), 8, 0, 1
  },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogY.Abs), 0, 0, 2 },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogY.Compare), 8, 0,
    1 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02XAngleWatchdog.Abs), 0, 0, 2 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02XAngleWatchdog.Compare), 8, 0, 1
  },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogX.Abs), 0, 0, 2 },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogX.Compare), 8, 0,
    1 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02YAngleWatchdog_l.Abs), 0, 0, 2 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02YAngleWatchdog_l.Compare), 8, 0,
    1 },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogY_a.Abs), 0, 0, 2
  },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogY_a.Compare), 8,
    0, 1 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02XAngleWatchdog_n.Abs), 0, 0, 2 },

  { (char_T *)(&ExperimentRobustControl_B.SRV02XAngleWatchdog_n.Compare), 8, 0,
    1 },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogX_d.Abs), 0, 0, 2
  },

  { (char_T *)(&ExperimentRobustControl_B.PendulumAngleWatchdogX_d.Compare), 8,
    0, 1 }
  ,

  { (char_T *)(&ExperimentRobustControl_DW.HILInitialize_AIMinimums[0]), 0, 0,
    64 },

  { (char_T *)(&ExperimentRobustControl_DW.HILInitialize_Card), 14, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_DW.HILReadEncoderTimebase_Task), 15, 0,
    1 },

  { (char_T *)(&ExperimentRobustControl_DW.HILWriteAnalog_PWORK), 11, 0, 7 },

  { (char_T *)(&ExperimentRobustControl_DW.HILInitialize_ClockModes[0]), 6, 0,
    59 },

  { (char_T *)(&ExperimentRobustControl_DW.HILInitialize_POSortedChans[0]), 7, 0,
    8 },

  { (char_T *)(&ExperimentRobustControl_DW.ResetEncoders_SubsysRanBC), 2, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_DW.HILInitialize_DOBits[0]), 8, 0, 8 },

  { (char_T *)(&ExperimentRobustControl_DW.HILWriteDigital_Buffer[0]), 16, 0, 4
  },

  { (char_T *)(&ExperimentRobustControl_DW.ResetEncoders_MODE), 8, 0, 1 }
};

/* data type transition table for block I/O structure */
static DataTypeTransitionTable rtBTransTable = {
  28U,
  rtBTransitions
};

/* data type transitions for Parameters structure */
static DataTypeTransition rtPTransitions[] = {
  { (char_T *)(&ExperimentRobustControl_P.ALPHA_MAX), 0, 0, 3280 },

  { (char_T *)(&ExperimentRobustControl_P.HILReadEncoderTimebase_clock), 6, 0,
    11 },

  { (char_T *)(&ExperimentRobustControl_P.HILInitialize_analog_input_chan[0]), 7,
    0, 56 },

  { (char_T *)(&ExperimentRobustControl_P.HILInitialize_active), 8, 0, 35 },

  { (char_T *)(&ExperimentRobustControl_P.StateSpace_X0), 0, 0, 69 },

  { (char_T *)(&ExperimentRobustControl_P.ManualSwitch2_CurrentSetting), 3, 0, 1
  },

  { (char_T *)(&ExperimentRobustControl_P.HILSetEncoderCounts_Active), 8, 0, 4 },

  { (char_T *)
    (&ExperimentRobustControl_P.SRV02YAngleWatchdog.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_P.SRV02YAngleWatchdog.Constant_Value), 0,
    0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.PendulumAngleWatchdogY.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_P.PendulumAngleWatchdogY.Constant_Value),
    0, 0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.SRV02XAngleWatchdog.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_P.SRV02XAngleWatchdog.Constant_Value), 0,
    0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.PendulumAngleWatchdogX.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_P.PendulumAngleWatchdogX.Constant_Value),
    0, 0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.SRV02YAngleWatchdog_l.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_P.SRV02YAngleWatchdog_l.Constant_Value),
    0, 0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.PendulumAngleWatchdogY_a.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.PendulumAngleWatchdogY_a.Constant_Value), 0, 0,
    1 },

  { (char_T *)
    (&ExperimentRobustControl_P.SRV02XAngleWatchdog_n.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)(&ExperimentRobustControl_P.SRV02XAngleWatchdog_n.Constant_Value),
    0, 0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.PendulumAngleWatchdogX_d.StopwithMessage_message_icon),
    2, 0, 1 },

  { (char_T *)
    (&ExperimentRobustControl_P.PendulumAngleWatchdogX_d.Constant_Value), 0, 0,
    1 }
};

/* data type transition table for Parameters structure */
static DataTypeTransitionTable rtPTransTable = {
  23U,
  rtPTransitions
};

/* [EOF] ExperimentRobustControl_dt.h */
