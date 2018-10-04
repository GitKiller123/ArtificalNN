/*
 * ExperimentRobustControl.h
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

#ifndef RTW_HEADER_ExperimentRobustControl_h_
#define RTW_HEADER_ExperimentRobustControl_h_
#include <math.h>
#include <string.h>
#ifndef ExperimentRobustControl_COMMON_INCLUDES_
# define ExperimentRobustControl_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "zero_crossing_types.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "dt_info.h"
#include "ext_work.h"
#include "hil.h"
#include "quanser_messages.h"
#include "quanser_extern.h"
#endif                                 /* ExperimentRobustControl_COMMON_INCLUDES_ */

#include "ExperimentRobustControl_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->blkStateChange = (val))
#endif

#ifndef rtmGetBlockIO
# define rtmGetBlockIO(rtm)            ((rtm)->blockIO)
#endif

#ifndef rtmSetBlockIO
# define rtmSetBlockIO(rtm, val)       ((rtm)->blockIO = (val))
#endif

#ifndef rtmGetChecksums
# define rtmGetChecksums(rtm)          ((rtm)->Sizes.checksums)
#endif

#ifndef rtmSetChecksums
# define rtmSetChecksums(rtm, val)     ((rtm)->Sizes.checksums = (val))
#endif

#ifndef rtmGetConstBlockIO
# define rtmGetConstBlockIO(rtm)       ((rtm)->constBlockIO)
#endif

#ifndef rtmSetConstBlockIO
# define rtmSetConstBlockIO(rtm, val)  ((rtm)->constBlockIO = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        ()
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)   ()
#endif

#ifndef rtmGetDefaultParam
# define rtmGetDefaultParam(rtm)       ((rtm)->defaultParam)
#endif

#ifndef rtmSetDefaultParam
# define rtmSetDefaultParam(rtm, val)  ((rtm)->defaultParam = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetDirectFeedThrough
# define rtmGetDirectFeedThrough(rtm)  ((rtm)->Sizes.sysDirFeedThru)
#endif

#ifndef rtmSetDirectFeedThrough
# define rtmSetDirectFeedThrough(rtm, val) ((rtm)->Sizes.sysDirFeedThru = (val))
#endif

#ifndef rtmGetErrorStatusFlag
# define rtmGetErrorStatusFlag(rtm)    ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatusFlag
# define rtmSetErrorStatusFlag(rtm, val) ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetFinalTime
# define rtmSetFinalTime(rtm, val)     ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetFirstInitCondFlag
# define rtmGetFirstInitCondFlag(rtm)  ()
#endif

#ifndef rtmSetFirstInitCondFlag
# define rtmSetFirstInitCondFlag(rtm, val) ()
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetMdlRefGlobalTID
# define rtmGetMdlRefGlobalTID(rtm)    ()
#endif

#ifndef rtmSetMdlRefGlobalTID
# define rtmSetMdlRefGlobalTID(rtm, val) ()
#endif

#ifndef rtmGetMdlRefTriggerTID
# define rtmGetMdlRefTriggerTID(rtm)   ()
#endif

#ifndef rtmSetMdlRefTriggerTID
# define rtmSetMdlRefTriggerTID(rtm, val) ()
#endif

#ifndef rtmGetModelMappingInfo
# define rtmGetModelMappingInfo(rtm)   ((rtm)->SpecialInfo.mappingInfo)
#endif

#ifndef rtmSetModelMappingInfo
# define rtmSetModelMappingInfo(rtm, val) ((rtm)->SpecialInfo.mappingInfo = (val))
#endif

#ifndef rtmGetModelName
# define rtmGetModelName(rtm)          ((rtm)->modelName)
#endif

#ifndef rtmSetModelName
# define rtmSetModelName(rtm, val)     ((rtm)->modelName = (val))
#endif

#ifndef rtmGetNonInlinedSFcns
# define rtmGetNonInlinedSFcns(rtm)    ()
#endif

#ifndef rtmSetNonInlinedSFcns
# define rtmSetNonInlinedSFcns(rtm, val) ()
#endif

#ifndef rtmGetNumBlockIO
# define rtmGetNumBlockIO(rtm)         ((rtm)->Sizes.numBlockIO)
#endif

#ifndef rtmSetNumBlockIO
# define rtmSetNumBlockIO(rtm, val)    ((rtm)->Sizes.numBlockIO = (val))
#endif

#ifndef rtmGetNumBlockParams
# define rtmGetNumBlockParams(rtm)     ((rtm)->Sizes.numBlockPrms)
#endif

#ifndef rtmSetNumBlockParams
# define rtmSetNumBlockParams(rtm, val) ((rtm)->Sizes.numBlockPrms = (val))
#endif

#ifndef rtmGetNumBlocks
# define rtmGetNumBlocks(rtm)          ((rtm)->Sizes.numBlocks)
#endif

#ifndef rtmSetNumBlocks
# define rtmSetNumBlocks(rtm, val)     ((rtm)->Sizes.numBlocks = (val))
#endif

#ifndef rtmGetNumContStates
# define rtmGetNumContStates(rtm)      ((rtm)->Sizes.numContStates)
#endif

#ifndef rtmSetNumContStates
# define rtmSetNumContStates(rtm, val) ((rtm)->Sizes.numContStates = (val))
#endif

#ifndef rtmGetNumDWork
# define rtmGetNumDWork(rtm)           ((rtm)->Sizes.numDwork)
#endif

#ifndef rtmSetNumDWork
# define rtmSetNumDWork(rtm, val)      ((rtm)->Sizes.numDwork = (val))
#endif

#ifndef rtmGetNumInputPorts
# define rtmGetNumInputPorts(rtm)      ((rtm)->Sizes.numIports)
#endif

#ifndef rtmSetNumInputPorts
# define rtmSetNumInputPorts(rtm, val) ((rtm)->Sizes.numIports = (val))
#endif

#ifndef rtmGetNumNonSampledZCs
# define rtmGetNumNonSampledZCs(rtm)   ((rtm)->Sizes.numNonSampZCs)
#endif

#ifndef rtmSetNumNonSampledZCs
# define rtmSetNumNonSampledZCs(rtm, val) ((rtm)->Sizes.numNonSampZCs = (val))
#endif

#ifndef rtmGetNumOutputPorts
# define rtmGetNumOutputPorts(rtm)     ((rtm)->Sizes.numOports)
#endif

#ifndef rtmSetNumOutputPorts
# define rtmSetNumOutputPorts(rtm, val) ((rtm)->Sizes.numOports = (val))
#endif

#ifndef rtmGetNumPeriodicContStates
# define rtmGetNumPeriodicContStates(rtm) ((rtm)->Sizes.numPeriodicContStates)
#endif

#ifndef rtmSetNumPeriodicContStates
# define rtmSetNumPeriodicContStates(rtm, val) ((rtm)->Sizes.numPeriodicContStates = (val))
#endif

#ifndef rtmGetNumSFcnParams
# define rtmGetNumSFcnParams(rtm)      ((rtm)->Sizes.numSFcnPrms)
#endif

#ifndef rtmSetNumSFcnParams
# define rtmSetNumSFcnParams(rtm, val) ((rtm)->Sizes.numSFcnPrms = (val))
#endif

#ifndef rtmGetNumSFunctions
# define rtmGetNumSFunctions(rtm)      ((rtm)->Sizes.numSFcns)
#endif

#ifndef rtmSetNumSFunctions
# define rtmSetNumSFunctions(rtm, val) ((rtm)->Sizes.numSFcns = (val))
#endif

#ifndef rtmGetNumSampleTimes
# define rtmGetNumSampleTimes(rtm)     ((rtm)->Sizes.numSampTimes)
#endif

#ifndef rtmSetNumSampleTimes
# define rtmSetNumSampleTimes(rtm, val) ((rtm)->Sizes.numSampTimes = (val))
#endif

#ifndef rtmGetNumU
# define rtmGetNumU(rtm)               ((rtm)->Sizes.numU)
#endif

#ifndef rtmSetNumU
# define rtmSetNumU(rtm, val)          ((rtm)->Sizes.numU = (val))
#endif

#ifndef rtmGetNumY
# define rtmGetNumY(rtm)               ((rtm)->Sizes.numY)
#endif

#ifndef rtmSetNumY
# define rtmSetNumY(rtm, val)          ((rtm)->Sizes.numY = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ()
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ()
#endif

#ifndef rtmGetOffsetTimeArray
# define rtmGetOffsetTimeArray(rtm)    ((rtm)->Timing.offsetTimesArray)
#endif

#ifndef rtmSetOffsetTimeArray
# define rtmSetOffsetTimeArray(rtm, val) ((rtm)->Timing.offsetTimesArray = (val))
#endif

#ifndef rtmGetOffsetTimePtr
# define rtmGetOffsetTimePtr(rtm)      ((rtm)->Timing.offsetTimes)
#endif

#ifndef rtmSetOffsetTimePtr
# define rtmSetOffsetTimePtr(rtm, val) ((rtm)->Timing.offsetTimes = (val))
#endif

#ifndef rtmGetOptions
# define rtmGetOptions(rtm)            ((rtm)->Sizes.options)
#endif

#ifndef rtmSetOptions
# define rtmSetOptions(rtm, val)       ((rtm)->Sizes.options = (val))
#endif

#ifndef rtmGetParamIsMalloced
# define rtmGetParamIsMalloced(rtm)    ()
#endif

#ifndef rtmSetParamIsMalloced
# define rtmSetParamIsMalloced(rtm, val) ()
#endif

#ifndef rtmGetPath
# define rtmGetPath(rtm)               ((rtm)->path)
#endif

#ifndef rtmSetPath
# define rtmSetPath(rtm, val)          ((rtm)->path = (val))
#endif

#ifndef rtmGetPerTaskSampleHits
# define rtmGetPerTaskSampleHits(rtm)  ()
#endif

#ifndef rtmSetPerTaskSampleHits
# define rtmSetPerTaskSampleHits(rtm, val) ()
#endif

#ifndef rtmGetPerTaskSampleHitsArray
# define rtmGetPerTaskSampleHitsArray(rtm) ((rtm)->Timing.perTaskSampleHitsArray)
#endif

#ifndef rtmSetPerTaskSampleHitsArray
# define rtmSetPerTaskSampleHitsArray(rtm, val) ((rtm)->Timing.perTaskSampleHitsArray = (val))
#endif

#ifndef rtmGetPerTaskSampleHitsPtr
# define rtmGetPerTaskSampleHitsPtr(rtm) ((rtm)->Timing.perTaskSampleHits)
#endif

#ifndef rtmSetPerTaskSampleHitsPtr
# define rtmSetPerTaskSampleHitsPtr(rtm, val) ((rtm)->Timing.perTaskSampleHits = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetPrevZCSigState
# define rtmGetPrevZCSigState(rtm)     ((rtm)->prevZCSigState)
#endif

#ifndef rtmSetPrevZCSigState
# define rtmSetPrevZCSigState(rtm, val) ((rtm)->prevZCSigState = (val))
#endif

#ifndef rtmGetRTWExtModeInfo
# define rtmGetRTWExtModeInfo(rtm)     ((rtm)->extModeInfo)
#endif

#ifndef rtmSetRTWExtModeInfo
# define rtmSetRTWExtModeInfo(rtm, val) ((rtm)->extModeInfo = (val))
#endif

#ifndef rtmGetRTWGeneratedSFcn
# define rtmGetRTWGeneratedSFcn(rtm)   ((rtm)->Sizes.rtwGenSfcn)
#endif

#ifndef rtmSetRTWGeneratedSFcn
# define rtmSetRTWGeneratedSFcn(rtm, val) ((rtm)->Sizes.rtwGenSfcn = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ()
#endif

#ifndef rtmSetRTWLogInfo
# define rtmSetRTWLogInfo(rtm, val)    ()
#endif

#ifndef rtmGetRTWRTModelMethodsInfo
# define rtmGetRTWRTModelMethodsInfo(rtm) ()
#endif

#ifndef rtmSetRTWRTModelMethodsInfo
# define rtmSetRTWRTModelMethodsInfo(rtm, val) ()
#endif

#ifndef rtmGetRTWSfcnInfo
# define rtmGetRTWSfcnInfo(rtm)        ((rtm)->sfcnInfo)
#endif

#ifndef rtmSetRTWSfcnInfo
# define rtmSetRTWSfcnInfo(rtm, val)   ((rtm)->sfcnInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfo
# define rtmGetRTWSolverInfo(rtm)      ((rtm)->solverInfo)
#endif

#ifndef rtmSetRTWSolverInfo
# define rtmSetRTWSolverInfo(rtm, val) ((rtm)->solverInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfoPtr
# define rtmGetRTWSolverInfoPtr(rtm)   ((rtm)->solverInfoPtr)
#endif

#ifndef rtmSetRTWSolverInfoPtr
# define rtmSetRTWSolverInfoPtr(rtm, val) ((rtm)->solverInfoPtr = (val))
#endif

#ifndef rtmGetReservedForXPC
# define rtmGetReservedForXPC(rtm)     ((rtm)->SpecialInfo.xpcData)
#endif

#ifndef rtmSetReservedForXPC
# define rtmSetReservedForXPC(rtm, val) ((rtm)->SpecialInfo.xpcData = (val))
#endif

#ifndef rtmGetRootDWork
# define rtmGetRootDWork(rtm)          ((rtm)->dwork)
#endif

#ifndef rtmSetRootDWork
# define rtmSetRootDWork(rtm, val)     ((rtm)->dwork = (val))
#endif

#ifndef rtmGetSFunctions
# define rtmGetSFunctions(rtm)         ((rtm)->childSfunctions)
#endif

#ifndef rtmSetSFunctions
# define rtmSetSFunctions(rtm, val)    ((rtm)->childSfunctions = (val))
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmSetSampleHitArray
# define rtmSetSampleHitArray(rtm, val) ((rtm)->Timing.sampleHitArray = (val))
#endif

#ifndef rtmGetSampleHitPtr
# define rtmGetSampleHitPtr(rtm)       ((rtm)->Timing.sampleHits)
#endif

#ifndef rtmSetSampleHitPtr
# define rtmSetSampleHitPtr(rtm, val)  ((rtm)->Timing.sampleHits = (val))
#endif

#ifndef rtmGetSampleTimeArray
# define rtmGetSampleTimeArray(rtm)    ((rtm)->Timing.sampleTimesArray)
#endif

#ifndef rtmSetSampleTimeArray
# define rtmSetSampleTimeArray(rtm, val) ((rtm)->Timing.sampleTimesArray = (val))
#endif

#ifndef rtmGetSampleTimePtr
# define rtmGetSampleTimePtr(rtm)      ((rtm)->Timing.sampleTimes)
#endif

#ifndef rtmSetSampleTimePtr
# define rtmSetSampleTimePtr(rtm, val) ((rtm)->Timing.sampleTimes = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDArray
# define rtmGetSampleTimeTaskIDArray(rtm) ((rtm)->Timing.sampleTimeTaskIDArray)
#endif

#ifndef rtmSetSampleTimeTaskIDArray
# define rtmSetSampleTimeTaskIDArray(rtm, val) ((rtm)->Timing.sampleTimeTaskIDArray = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDPtr
# define rtmGetSampleTimeTaskIDPtr(rtm) ((rtm)->Timing.sampleTimeTaskIDPtr)
#endif

#ifndef rtmSetSampleTimeTaskIDPtr
# define rtmSetSampleTimeTaskIDPtr(rtm, val) ((rtm)->Timing.sampleTimeTaskIDPtr = (val))
#endif

#ifndef rtmGetSelf
# define rtmGetSelf(rtm)               ()
#endif

#ifndef rtmSetSelf
# define rtmSetSelf(rtm, val)          ()
#endif

#ifndef rtmGetSimMode
# define rtmGetSimMode(rtm)            ((rtm)->simMode)
#endif

#ifndef rtmSetSimMode
# define rtmSetSimMode(rtm, val)       ((rtm)->simMode = (val))
#endif

#ifndef rtmGetSimTimeStep
# define rtmGetSimTimeStep(rtm)        ((rtm)->Timing.simTimeStep)
#endif

#ifndef rtmSetSimTimeStep
# define rtmSetSimTimeStep(rtm, val)   ((rtm)->Timing.simTimeStep = (val))
#endif

#ifndef rtmGetStartTime
# define rtmGetStartTime(rtm)          ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetStartTime
# define rtmSetStartTime(rtm, val)     ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmSetStepSize
# define rtmSetStepSize(rtm, val)      ((rtm)->Timing.stepSize = (val))
#endif

#ifndef rtmGetStopRequestedFlag
# define rtmGetStopRequestedFlag(rtm)  ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequestedFlag
# define rtmSetStopRequestedFlag(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetTaskCounters
# define rtmGetTaskCounters(rtm)       ()
#endif

#ifndef rtmSetTaskCounters
# define rtmSetTaskCounters(rtm, val)  ()
#endif

#ifndef rtmGetTaskTimeArray
# define rtmGetTaskTimeArray(rtm)      ((rtm)->Timing.tArray)
#endif

#ifndef rtmSetTaskTimeArray
# define rtmSetTaskTimeArray(rtm, val) ((rtm)->Timing.tArray = (val))
#endif

#ifndef rtmGetTimePtr
# define rtmGetTimePtr(rtm)            ((rtm)->Timing.t)
#endif

#ifndef rtmSetTimePtr
# define rtmSetTimePtr(rtm, val)       ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTimingData
# define rtmGetTimingData(rtm)         ((rtm)->Timing.timingData)
#endif

#ifndef rtmSetTimingData
# define rtmSetTimingData(rtm, val)    ((rtm)->Timing.timingData = (val))
#endif

#ifndef rtmGetU
# define rtmGetU(rtm)                  ((rtm)->inputs)
#endif

#ifndef rtmSetU
# define rtmSetU(rtm, val)             ((rtm)->inputs = (val))
#endif

#ifndef rtmGetVarNextHitTimesListPtr
# define rtmGetVarNextHitTimesListPtr(rtm) ((rtm)->Timing.varNextHitTimesList)
#endif

#ifndef rtmSetVarNextHitTimesListPtr
# define rtmSetVarNextHitTimesListPtr(rtm, val) ((rtm)->Timing.varNextHitTimesList = (val))
#endif

#ifndef rtmGetY
# define rtmGetY(rtm)                  ((rtm)->outputs)
#endif

#ifndef rtmSetY
# define rtmSetY(rtm, val)             ((rtm)->outputs = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetZCSignalValues
# define rtmGetZCSignalValues(rtm)     ((rtm)->zcSignalValues)
#endif

#ifndef rtmSetZCSignalValues
# define rtmSetZCSignalValues(rtm, val) ((rtm)->zcSignalValues = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmSet_TimeOfLastOutput
# define rtmSet_TimeOfLastOutput(rtm, val) ((rtm)->Timing.timeOfLastOutput = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetChecksumVal
# define rtmGetChecksumVal(rtm, idx)   ((rtm)->Sizes.checksums[idx])
#endif

#ifndef rtmSetChecksumVal
# define rtmSetChecksumVal(rtm, idx, val) ((rtm)->Sizes.checksums[idx] = (val))
#endif

#ifndef rtmGetDWork
# define rtmGetDWork(rtm, idx)         ((rtm)->dwork[idx])
#endif

#ifndef rtmSetDWork
# define rtmSetDWork(rtm, idx, val)    ((rtm)->dwork[idx] = (val))
#endif

#ifndef rtmGetOffsetTime
# define rtmGetOffsetTime(rtm, idx)    ((rtm)->Timing.offsetTimes[idx])
#endif

#ifndef rtmSetOffsetTime
# define rtmSetOffsetTime(rtm, idx, val) ((rtm)->Timing.offsetTimes[idx] = (val))
#endif

#ifndef rtmGetSFunction
# define rtmGetSFunction(rtm, idx)     ((rtm)->childSfunctions[idx])
#endif

#ifndef rtmSetSFunction
# define rtmSetSFunction(rtm, idx, val) ((rtm)->childSfunctions[idx] = (val))
#endif

#ifndef rtmGetSampleTime
# define rtmGetSampleTime(rtm, idx)    ((rtm)->Timing.sampleTimes[idx])
#endif

#ifndef rtmSetSampleTime
# define rtmSetSampleTime(rtm, idx, val) ((rtm)->Timing.sampleTimes[idx] = (val))
#endif

#ifndef rtmGetSampleTimeTaskID
# define rtmGetSampleTimeTaskID(rtm, idx) ((rtm)->Timing.sampleTimeTaskIDPtr[idx])
#endif

#ifndef rtmSetSampleTimeTaskID
# define rtmSetSampleTimeTaskID(rtm, idx, val) ((rtm)->Timing.sampleTimeTaskIDPtr[idx] = (val))
#endif

#ifndef rtmGetVarNextHitTimeList
# define rtmGetVarNextHitTimeList(rtm, idx) ((rtm)->Timing.varNextHitTimesList[idx])
#endif

#ifndef rtmSetVarNextHitTimeList
# define rtmSetVarNextHitTimeList(rtm, idx, val) ((rtm)->Timing.varNextHitTimesList[idx] = (val))
#endif

#ifndef rtmIsContinuousTask
# define rtmIsContinuousTask(rtm, tid) ((tid) == 0)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmIsSampleHit
# define rtmIsSampleHit(rtm, sti, tid) ((rtmIsMajorTimeStep((rtm)) && (rtm)->Timing.sampleHits[(rtm)->Timing.sampleTimeTaskIDPtr[sti]]))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmSetT
# define rtmSetT(rtm, val)                                       /* Do Nothing */
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetTStart
# define rtmSetTStart(rtm, val)        ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetTaskTime
# define rtmGetTaskTime(rtm, sti)      (rtmGetTPtr((rtm))[(rtm)->Timing.sampleTimeTaskIDPtr[sti]])
#endif

#ifndef rtmSetTaskTime
# define rtmSetTaskTime(rtm, sti, val) (rtmGetTPtr((rtm))[sti] = (val))
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifdef rtmGetRTWSolverInfo
#undef rtmGetRTWSolverInfo
#endif

#define rtmGetRTWSolverInfo(rtm)       &((rtm)->solverInfo)

/* Definition for use in the target main file */
#define ExperimentRobustControl_rtModel RT_MODEL_ExperimentRobustCont_T

/* Block signals for system '<S6>/Pendulum Angle Watchdog: X' */
typedef struct {
  real_T Abs;                          /* '<S9>/Abs' */
  real_T CheckMaximumAllowedAngle;     /* '<S9>/Check Maximum Allowed Angle' */
  boolean_T Compare;                   /* '<S12>/Compare' */
} B_PendulumAngleWatchdogX_Expe_T;

/* Block signals for system '<S6>/SRV02 X Angle Watchdog' */
typedef struct {
  real_T Abs;                          /* '<S10>/Abs' */
  real_T CheckMaximumAllowedAngle;     /* '<S10>/Check Maximum Allowed Angle' */
  boolean_T Compare;                   /* '<S14>/Compare' */
} B_SRV02XAngleWatchdog_Experim_T;

/* Block signals for system '<S7>/Pendulum Angle Watchdog: Y' */
typedef struct {
  real_T Abs;                          /* '<S15>/Abs' */
  real_T CheckMaximumAllowedAngle;     /* '<S15>/Check Maximum Allowed Angle' */
  boolean_T Compare;                   /* '<S18>/Compare' */
} B_PendulumAngleWatchdogY_Expe_T;

/* Block signals for system '<S7>/SRV02 Y Angle Watchdog' */
typedef struct {
  real_T Abs;                          /* '<S16>/Abs' */
  real_T CheckMaximumAllowedAngle;     /* '<S16>/Check Maximum Allowed Angle' */
  boolean_T Compare;                   /* '<S20>/Compare' */
} B_SRV02YAngleWatchdog_Experim_T;

/* Block signals (auto storage) */
typedef struct {
  real_T HILReadEncoderTimebase_o1;    /* '<S3>/HIL Read Encoder Timebase' */
  real_T HILReadEncoderTimebase_o2;    /* '<S3>/HIL Read Encoder Timebase' */
  real_T HILReadEncoderTimebase_o3;    /* '<S3>/HIL Read Encoder Timebase' */
  real_T HILReadEncoderTimebase_o4;    /* '<S3>/HIL Read Encoder Timebase' */
  real_T uDIPEncoderCalibrationradcount;/* '<S23>/2D IP Encoder Calibration (rad//count)' */
  real_T EncoderCalibrationradcount;   /* '<S23>/Encoder Calibration  (rad//count)' */
  real_T uDIPEncoderYCalibrationradcount;/* '<S24>/2D IP Encoder Y  Calibration (rad//count)' */
  real_T EncoderCalibrationYradcount;  /* '<S24>/Encoder Calibration Y (rad//count)' */
  real_T StateSpace[2];                /* '<S1>/State-Space' */
  real_T Integrator[4];                /* '<S4>/Integrator' */
  real_T Integrator_m;                 /* '<S6>/Integrator' */
  real_T Integrator_k[4];              /* '<S5>/Integrator' */
  real_T Integrator_c;                 /* '<S7>/Integrator' */
  real_T TmpSignalConversionAtLQRGainYIn[5];
  real_T ManualSwitch2[2];             /* '<Root>/Manual Switch2' */
  real_T Gain1[2];                     /* '<Root>/Gain1' */
  real_T DirectionConventionRightHandsys[2];/* '<S21>/Direction Convention: (Right-Hand) system' */
  real_T AmplifierSaturationV[2];      /* '<S21>/Amplifier  Saturation (V)' */
  real_T InverseAmplifierGainVV[2];    /* '<S21>/Inverse Amplifier  Gain (V//V)' */
  real_T DACBSaturationV[2];           /* '<S21>/DACB Saturation (V)' */
  real_T AmplifierGainVV[2];           /* '<S21>/Amplifier  Gain (V//V)' */
  real_T Clock;                        /* '<S3>/Clock' */
  real_T Gain2[2];                     /* '<Root>/Gain2' */
  real_T Gain[2];                      /* '<Root>/Gain' */
  real_T TmpSignalConversionAtStateSpace[4];
  real_T Gain_l[5];                    /* '<S8>/Gain' */
  real_T syscc_x[2];                   /* '<S4>/-sysc.c_x' */
  real_T Sum1[2];                      /* '<S4>/Sum1' */
  real_T KalmanFilterGain[4];          /* '<S4>/Kalman Filter Gain' */
  real_T syscb_x[4];                   /* '<S4>/sysc.b_x' */
  real_T sysca_x[4];                   /* '<S4>/sysc.a_x' */
  real_T Sum[4];                       /* '<S4>/Sum' */
  real_T sysdc_x[2];                   /* '<S5>/-sysd.c_x' */
  real_T Sum1_k[2];                    /* '<S5>/Sum1' */
  real_T KalmanFilterGain_e[4];        /* '<S5>/Kalman Filter Gain' */
  real_T sysdb_x[4];                   /* '<S5>/sysd.b_x' */
  real_T sysda_x[4];                   /* '<S5>/sysd.a_x' */
  real_T Sum_j[4];                     /* '<S5>/Sum' */
  real_T TmpSignalConversionAtLQRGainXIn[5];
  real_T LQRGainX;                     /* '<S2>/LQR Gain: X' */
  real_T LQRGainY;                     /* '<S2>/LQR Gain: Y' */
  uint8_T Compare;                     /* '<S22>/Compare' */
  B_SRV02YAngleWatchdog_Experim_T SRV02YAngleWatchdog;/* '<S24>/SRV02 Y Angle Watchdog' */
  B_PendulumAngleWatchdogY_Expe_T PendulumAngleWatchdogY;/* '<S24>/Pendulum Angle Watchdog: Y' */
  B_SRV02XAngleWatchdog_Experim_T SRV02XAngleWatchdog;/* '<S23>/SRV02 X Angle Watchdog' */
  B_PendulumAngleWatchdogX_Expe_T PendulumAngleWatchdogX;/* '<S23>/Pendulum Angle Watchdog: X' */
  B_SRV02YAngleWatchdog_Experim_T SRV02YAngleWatchdog_l;/* '<S7>/SRV02 Y Angle Watchdog' */
  B_PendulumAngleWatchdogY_Expe_T PendulumAngleWatchdogY_a;/* '<S7>/Pendulum Angle Watchdog: Y' */
  B_SRV02XAngleWatchdog_Experim_T SRV02XAngleWatchdog_n;/* '<S6>/SRV02 X Angle Watchdog' */
  B_PendulumAngleWatchdogX_Expe_T PendulumAngleWatchdogX_d;/* '<S6>/Pendulum Angle Watchdog: X' */
} B_ExperimentRobustControl_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T HILInitialize_AIMinimums[8];  /* '<S3>/HIL Initialize' */
  real_T HILInitialize_AIMaximums[8];  /* '<S3>/HIL Initialize' */
  real_T HILInitialize_AOMinimums[8];  /* '<S3>/HIL Initialize' */
  real_T HILInitialize_AOMaximums[8];  /* '<S3>/HIL Initialize' */
  real_T HILInitialize_AOVoltages[8];  /* '<S3>/HIL Initialize' */
  real_T HILInitialize_FilterFrequency[8];/* '<S3>/HIL Initialize' */
  real_T HILInitialize_POSortedFreqs[8];/* '<S3>/HIL Initialize' */
  real_T HILInitialize_POValues[8];    /* '<S3>/HIL Initialize' */
  t_card HILInitialize_Card;           /* '<S3>/HIL Initialize' */
  t_task HILReadEncoderTimebase_Task;  /* '<S3>/HIL Read Encoder Timebase' */
  void *HILWriteAnalog_PWORK;          /* '<S3>/HIL Write Analog' */
  void *HILWriteDigital_PWORK;         /* '<S3>/HIL Write Digital' */
  struct {
    void *LoggedData[3];
  } Scope1_PWORK;                      /* '<Root>/Scope1' */

  struct {
    void *LoggedData;
  } Scope_PWORK;                       /* '<S2>/Scope' */

  void *HILSetEncoderCounts_PWORK;     /* '<S25>/HIL Set Encoder Counts' */
  int32_T HILInitialize_ClockModes[3]; /* '<S3>/HIL Initialize' */
  int32_T HILInitialize_DOStates[8];   /* '<S3>/HIL Initialize' */
  int32_T HILInitialize_QuadratureModes[8];/* '<S3>/HIL Initialize' */
  int32_T HILInitialize_InitialEICounts[8];/* '<S3>/HIL Initialize' */
  int32_T HILInitialize_POModeValues[8];/* '<S3>/HIL Initialize' */
  int32_T HILInitialize_POAlignValues[8];/* '<S3>/HIL Initialize' */
  int32_T HILInitialize_POPolarityVals[8];/* '<S3>/HIL Initialize' */
  int32_T HILReadEncoderTimebase_Buffer[4];/* '<S3>/HIL Read Encoder Timebase' */
  int32_T HILSetEncoderCounts_Buffer[4];/* '<S25>/HIL Set Encoder Counts' */
  uint32_T HILInitialize_POSortedChans[8];/* '<S3>/HIL Initialize' */
  int8_T ResetEncoders_SubsysRanBC;    /* '<S3>/Reset Encoders' */
  boolean_T HILInitialize_DOBits[8];   /* '<S3>/HIL Initialize' */
  t_boolean HILWriteDigital_Buffer[4]; /* '<S3>/HIL Write Digital' */
  boolean_T ResetEncoders_MODE;        /* '<S3>/Reset Encoders' */
} DW_ExperimentRobustControl_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T StateSpace_CSTATE[54];        /* '<S1>/State-Space' */
  real_T Integrator_CSTATE[4];         /* '<S4>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S6>/Integrator' */
  real_T Integrator_CSTATE_n[4];       /* '<S5>/Integrator' */
  real_T Integrator_CSTATE_a;          /* '<S7>/Integrator' */
} X_ExperimentRobustControl_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T StateSpace_CSTATE[54];        /* '<S1>/State-Space' */
  real_T Integrator_CSTATE[4];         /* '<S4>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S6>/Integrator' */
  real_T Integrator_CSTATE_n[4];       /* '<S5>/Integrator' */
  real_T Integrator_CSTATE_a;          /* '<S7>/Integrator' */
} XDot_ExperimentRobustControl_T;

/* State disabled  */
typedef struct {
  boolean_T StateSpace_CSTATE[54];     /* '<S1>/State-Space' */
  boolean_T Integrator_CSTATE[4];      /* '<S4>/Integrator' */
  boolean_T Integrator_CSTATE_c;       /* '<S6>/Integrator' */
  boolean_T Integrator_CSTATE_n[4];    /* '<S5>/Integrator' */
  boolean_T Integrator_CSTATE_a;       /* '<S7>/Integrator' */
} XDis_ExperimentRobustControl_T;

#ifndef ODE1_INTG
#define ODE1_INTG

/* ODE1 Integration Data */
typedef struct {
  real_T *f[1];                        /* derivatives */
} ODE1_IntgData;

#endif

/* Backward compatible GRT Identifiers */
#define rtB                            ExperimentRobustControl_B
#define BlockIO                        B_ExperimentRobustControl_T
#define rtX                            ExperimentRobustControl_X
#define ContinuousStates               X_ExperimentRobustControl_T
#define rtXdot                         ExperimentRobustControl_XDot
#define StateDerivatives               XDot_ExperimentRobustControl_T
#define tXdis                          ExperimentRobustControl_XDis
#define StateDisabled                  XDis_ExperimentRobustControl_T
#define rtP                            ExperimentRobustControl_P
#define Parameters                     P_ExperimentRobustControl_T
#define rtDWork                        ExperimentRobustControl_DW
#define D_Work                         DW_ExperimentRobustControl_T

/* Parameters for system: '<S6>/Pendulum Angle Watchdog: X' */
struct P_PendulumAngleWatchdogX_Expe_T_ {
  int8_T StopwithMessage_message_icon; /* Mask Parameter: StopwithMessage_message_icon
                                        * Referenced by: '<S11>/Show Message on Host'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S12>/Constant'
                                        */
};

/* Parameters for system: '<S6>/SRV02 X Angle Watchdog' */
struct P_SRV02XAngleWatchdog_Experim_T_ {
  int8_T StopwithMessage_message_icon; /* Mask Parameter: StopwithMessage_message_icon
                                        * Referenced by: '<S13>/Show Message on Host'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S14>/Constant'
                                        */
};

/* Parameters for system: '<S7>/Pendulum Angle Watchdog: Y' */
struct P_PendulumAngleWatchdogY_Expe_T_ {
  int8_T StopwithMessage_message_icon; /* Mask Parameter: StopwithMessage_message_icon
                                        * Referenced by: '<S17>/Show Message on Host'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S18>/Constant'
                                        */
};

/* Parameters for system: '<S7>/SRV02 Y Angle Watchdog' */
struct P_SRV02YAngleWatchdog_Experim_T_ {
  int8_T StopwithMessage_message_icon; /* Mask Parameter: StopwithMessage_message_icon
                                        * Referenced by: '<S19>/Show Message on Host'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S20>/Constant'
                                        */
};

/* Parameters (auto storage) */
struct P_ExperimentRobustControl_T_ {
  real_T ALPHA_MAX;                    /* Variable: ALPHA_MAX
                                        * Referenced by:
                                        *   '<S9>/X MAX'
                                        *   '<S15>/X MAX'
                                        *   '<S26>/X MAX'
                                        *   '<S32>/X MAX'
                                        */
  real_T K_AMP;                        /* Variable: K_AMP
                                        * Referenced by:
                                        *   '<S21>/Amplifier  Gain (V//V)'
                                        *   '<S21>/Inverse Amplifier  Gain (V//V)'
                                        */
  real_T K_CT_labb[5];                 /* Variable: K_CT_labb
                                        * Referenced by:
                                        *   '<S2>/LQR Gain: X'
                                        *   '<S2>/LQR Gain: Y'
                                        */
  real_T K_ENC;                        /* Variable: K_ENC
                                        * Referenced by:
                                        *   '<S23>/Encoder Calibration  (rad//count)'
                                        *   '<S24>/Encoder Calibration Y (rad//count)'
                                        */
  real_T K_ENC_2DIP[2];                /* Variable: K_ENC_2DIP
                                        * Referenced by:
                                        *   '<S23>/2D IP Encoder Calibration (rad//count)'
                                        *   '<S24>/2D IP Encoder Y  Calibration (rad//count)'
                                        */
  real_T L_CT_labb[8];                 /* Variable: L_CT_labb
                                        * Referenced by:
                                        *   '<S4>/Kalman Filter Gain'
                                        *   '<S5>/Kalman Filter Gain'
                                        */
  real_T OCL;                          /* Variable: OCL
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T THETA_MAX;                    /* Variable: THETA_MAX
                                        * Referenced by:
                                        *   '<S10>/X MAX'
                                        *   '<S16>/X MAX'
                                        *   '<S27>/X MAX'
                                        *   '<S33>/X MAX'
                                        */
  real_T VMAX_AMP;                     /* Variable: VMAX_AMP
                                        * Referenced by: '<S21>/Amplifier  Saturation (V)'
                                        */
  real_T VMAX_DAC;                     /* Variable: VMAX_DAC
                                        * Referenced by: '<S21>/DACB Saturation (V)'
                                        */
  real_T ah[2916];                     /* Variable: ah
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T bh[216];                      /* Variable: bh
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T ch[108];                      /* Variable: ch
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T HILInitialize_analog_input_maxi;/* Mask Parameter: HILInitialize_analog_input_maxi
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_analog_input_mini;/* Mask Parameter: HILInitialize_analog_input_mini
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_analog_output_max;/* Mask Parameter: HILInitialize_analog_output_max
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_analog_output_min;/* Mask Parameter: HILInitialize_analog_output_min
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S22>/Constant'
                                        */
  real_T HILInitialize_final_analog_outp;/* Mask Parameter: HILInitialize_final_analog_outp
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_final_pwm_outputs;/* Mask Parameter: HILInitialize_final_pwm_outputs
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_initial_analog_ou;/* Mask Parameter: HILInitialize_initial_analog_ou
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_initial_pwm_outpu;/* Mask Parameter: HILInitialize_initial_pwm_outpu
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_pwm_frequency;  /* Mask Parameter: HILInitialize_pwm_frequency
                                        * Referenced by: '<S3>/HIL Initialize'
                                        */
  real_T HILInitialize_pwm_leading_deadb;/* Mask Parameter: HILInitialize_pwm_leading_deadb
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_pwm_trailing_dead;/* Mask Parameter: HILInitialize_pwm_trailing_dead
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_set_other_outputs;/* Mask Parameter: HILInitialize_set_other_outputs
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_set_other_outpu_c;/* Mask Parameter: HILInitialize_set_other_outpu_c
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_set_other_outpu_p;/* Mask Parameter: HILInitialize_set_other_outpu_p
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_set_other_outpu_h;/* Mask Parameter: HILInitialize_set_other_outpu_h
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_watchdog_analog_o;/* Mask Parameter: HILInitialize_watchdog_analog_o
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  real_T HILInitialize_watchdog_pwm_outp;/* Mask Parameter: HILInitialize_watchdog_pwm_outp
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  int32_T HILReadEncoderTimebase_clock;/* Mask Parameter: HILReadEncoderTimebase_clock
                                        * Referenced by: '<S3>/HIL Read Encoder Timebase'
                                        */
  int32_T HILSetEncoderCounts_counts;  /* Mask Parameter: HILSetEncoderCounts_counts
                                        * Referenced by: '<S25>/HIL Set Encoder Counts'
                                        */
  int32_T HILInitialize_hardware_clocks[3];/* Mask Parameter: HILInitialize_hardware_clocks
                                            * Referenced by: '<S3>/HIL Initialize'
                                            */
  int32_T HILInitialize_initial_encoder_c;/* Mask Parameter: HILInitialize_initial_encoder_c
                                           * Referenced by: '<S3>/HIL Initialize'
                                           */
  int32_T HILInitialize_pwm_alignment; /* Mask Parameter: HILInitialize_pwm_alignment
                                        * Referenced by: '<S3>/HIL Initialize'
                                        */
  int32_T HILInitialize_pwm_configuration;/* Mask Parameter: HILInitialize_pwm_configuration
                                           * Referenced by: '<S3>/HIL Initialize'
                                           */
  int32_T HILInitialize_pwm_modes;     /* Mask Parameter: HILInitialize_pwm_modes
                                        * Referenced by: '<S3>/HIL Initialize'
                                        */
  int32_T HILInitialize_pwm_polarity;  /* Mask Parameter: HILInitialize_pwm_polarity
                                        * Referenced by: '<S3>/HIL Initialize'
                                        */
  int32_T HILInitialize_watchdog_digital_;/* Mask Parameter: HILInitialize_watchdog_digital_
                                           * Referenced by: '<S3>/HIL Initialize'
                                           */
  uint32_T HILInitialize_analog_input_chan[8];/* Mask Parameter: HILInitialize_analog_input_chan
                                               * Referenced by: '<S3>/HIL Initialize'
                                               */
  uint32_T HILInitialize_analog_output_cha[8];/* Mask Parameter: HILInitialize_analog_output_cha
                                               * Referenced by: '<S3>/HIL Initialize'
                                               */
  uint32_T HILSetEncoderCounts_channels[4];/* Mask Parameter: HILSetEncoderCounts_channels
                                            * Referenced by: '<S25>/HIL Set Encoder Counts'
                                            */
  uint32_T HILReadEncoderTimebase_channels[4];/* Mask Parameter: HILReadEncoderTimebase_channels
                                               * Referenced by: '<S3>/HIL Read Encoder Timebase'
                                               */
  uint32_T HILWriteAnalog_channels[2]; /* Mask Parameter: HILWriteAnalog_channels
                                        * Referenced by: '<S3>/HIL Write Analog'
                                        */
  uint32_T HILWriteDigital_channels[4];/* Mask Parameter: HILWriteDigital_channels
                                        * Referenced by: '<S3>/HIL Write Digital'
                                        */
  uint32_T HILInitialize_digital_output_ch[8];/* Mask Parameter: HILInitialize_digital_output_ch
                                               * Referenced by: '<S3>/HIL Initialize'
                                               */
  uint32_T HILInitialize_encoder_channels[8];/* Mask Parameter: HILInitialize_encoder_channels
                                              * Referenced by: '<S3>/HIL Initialize'
                                              */
  uint32_T HILInitialize_pwm_channels[8];/* Mask Parameter: HILInitialize_pwm_channels
                                          * Referenced by: '<S3>/HIL Initialize'
                                          */
  uint32_T HILInitialize_quadrature;   /* Mask Parameter: HILInitialize_quadrature
                                        * Referenced by: '<S3>/HIL Initialize'
                                        */
  uint32_T HILReadEncoderTimebase_samples_;/* Mask Parameter: HILReadEncoderTimebase_samples_
                                            * Referenced by: '<S3>/HIL Read Encoder Timebase'
                                            */
  boolean_T HILInitialize_active;      /* Mask Parameter: HILInitialize_active
                                        * Referenced by: '<S3>/HIL Initialize'
                                        */
  boolean_T HILInitialize_final_digital_out;/* Mask Parameter: HILInitialize_final_digital_out
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_initial_digital_o;/* Mask Parameter: HILInitialize_initial_digital_o
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_input_;/* Mask Parameter: HILInitialize_set_analog_input_
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_inpu_o;/* Mask Parameter: HILInitialize_set_analog_inpu_o
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_output;/* Mask Parameter: HILInitialize_set_analog_output
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_outp_d;/* Mask Parameter: HILInitialize_set_analog_outp_d
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_outp_j;/* Mask Parameter: HILInitialize_set_analog_outp_j
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_outp_m;/* Mask Parameter: HILInitialize_set_analog_outp_m
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_out_mi;/* Mask Parameter: HILInitialize_set_analog_out_mi
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_outp_n;/* Mask Parameter: HILInitialize_set_analog_outp_n
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_analog_outp_k;/* Mask Parameter: HILInitialize_set_analog_outp_k
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_clock_frequen;/* Mask Parameter: HILInitialize_set_clock_frequen
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_clock_frequ_l;/* Mask Parameter: HILInitialize_set_clock_frequ_l
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_clock_params_;/* Mask Parameter: HILInitialize_set_clock_params_
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_clock_param_o;/* Mask Parameter: HILInitialize_set_clock_param_o
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_digital_outpu;/* Mask Parameter: HILInitialize_set_digital_outpu
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_digital_out_j;/* Mask Parameter: HILInitialize_set_digital_out_j
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_digital_out_a;/* Mask Parameter: HILInitialize_set_digital_out_a
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_digital_out_f;/* Mask Parameter: HILInitialize_set_digital_out_f
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_digital_out_b;/* Mask Parameter: HILInitialize_set_digital_out_b
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_digital_ou_bw;/* Mask Parameter: HILInitialize_set_digital_ou_bw
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_digital_ou_j0;/* Mask Parameter: HILInitialize_set_digital_ou_j0
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_encoder_count;/* Mask Parameter: HILInitialize_set_encoder_count
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_encoder_cou_e;/* Mask Parameter: HILInitialize_set_encoder_cou_e
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_encoder_param;/* Mask Parameter: HILInitialize_set_encoder_param
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_encoder_par_m;/* Mask Parameter: HILInitialize_set_encoder_par_m
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_other_outpu_i;/* Mask Parameter: HILInitialize_set_other_outpu_i
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_pwm_outputs_a;/* Mask Parameter: HILInitialize_set_pwm_outputs_a
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_pwm_outputs_l;/* Mask Parameter: HILInitialize_set_pwm_outputs_l
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_pwm_outputs_n;/* Mask Parameter: HILInitialize_set_pwm_outputs_n
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_pwm_outputs_c;/* Mask Parameter: HILInitialize_set_pwm_outputs_c
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_pwm_outputs_o;/* Mask Parameter: HILInitialize_set_pwm_outputs_o
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_pwm_params_at;/* Mask Parameter: HILInitialize_set_pwm_params_at
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  boolean_T HILInitialize_set_pwm_params__m;/* Mask Parameter: HILInitialize_set_pwm_params__m
                                             * Referenced by: '<S3>/HIL Initialize'
                                             */
  real_T StateSpace_X0;                /* Expression: 0
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S4>/Integrator'
                                        */
  real_T Integrator_IC_b;              /* Expression: 0
                                        * Referenced by: '<S6>/Integrator'
                                        */
  real_T Integrator_IC_p;              /* Expression: 0
                                        * Referenced by: '<S5>/Integrator'
                                        */
  real_T Integrator_IC_px;             /* Expression: 0
                                        * Referenced by: '<S7>/Integrator'
                                        */
  real_T DirectionConventionRightHandsys;/* Expression: -1
                                          * Referenced by: '<S21>/Direction Convention: (Right-Hand) system'
                                          */
  real_T EnableVoltPAQX2X4_Value[4];   /* Expression: [1 1 1 1]
                                        * Referenced by: '<S3>/Enable VoltPAQ-X2,X4'
                                        */
  real_T Gain2_Gain;                   /* Expression: 180/pi
                                        * Referenced by: '<Root>/Gain2'
                                        */
  real_T Gain_Gain;                    /* Expression: 180/pi
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Gain_Gain_d;                  /* Expression: 180/pi
                                        * Referenced by: '<S8>/Gain'
                                        */
  real_T syscc_x_Gain[8];              /* Expression: sysc_labb.c
                                        * Referenced by: '<S4>/-sysc.c_x'
                                        */
  real_T syscb_x_Gain[4];              /* Expression: sysc_labb.b
                                        * Referenced by: '<S4>/sysc.b_x'
                                        */
  real_T sysca_x_Gain[16];             /* Expression: sysc_labb.a
                                        * Referenced by: '<S4>/sysc.a_x'
                                        */
  real_T sysdc_x_Gain[8];              /* Expression: sysc_labb.c
                                        * Referenced by: '<S5>/-sysd.c_x'
                                        */
  real_T sysdb_x_Gain[4];              /* Expression: sysc_labb.b
                                        * Referenced by: '<S5>/sysd.b_x'
                                        */
  real_T sysda_x_Gain[16];             /* Expression: sysc_labb.a
                                        * Referenced by: '<S5>/sysd.a_x'
                                        */
  uint8_T ManualSwitch2_CurrentSetting;/* Computed Parameter: ManualSwitch2_CurrentSetting
                                        * Referenced by: '<Root>/Manual Switch2'
                                        */
  boolean_T HILSetEncoderCounts_Active;/* Computed Parameter: HILSetEncoderCounts_Active
                                        * Referenced by: '<S25>/HIL Set Encoder Counts'
                                        */
  boolean_T HILReadEncoderTimebase_Active;/* Computed Parameter: HILReadEncoderTimebase_Active
                                           * Referenced by: '<S3>/HIL Read Encoder Timebase'
                                           */
  boolean_T HILWriteAnalog_Active;     /* Computed Parameter: HILWriteAnalog_Active
                                        * Referenced by: '<S3>/HIL Write Analog'
                                        */
  boolean_T HILWriteDigital_Active;    /* Computed Parameter: HILWriteDigital_Active
                                        * Referenced by: '<S3>/HIL Write Digital'
                                        */
  P_SRV02YAngleWatchdog_Experim_T SRV02YAngleWatchdog;/* '<S24>/SRV02 Y Angle Watchdog' */
  P_PendulumAngleWatchdogY_Expe_T PendulumAngleWatchdogY;/* '<S24>/Pendulum Angle Watchdog: Y' */
  P_SRV02XAngleWatchdog_Experim_T SRV02XAngleWatchdog;/* '<S23>/SRV02 X Angle Watchdog' */
  P_PendulumAngleWatchdogX_Expe_T PendulumAngleWatchdogX;/* '<S23>/Pendulum Angle Watchdog: X' */
  P_SRV02YAngleWatchdog_Experim_T SRV02YAngleWatchdog_l;/* '<S7>/SRV02 Y Angle Watchdog' */
  P_PendulumAngleWatchdogY_Expe_T PendulumAngleWatchdogY_a;/* '<S7>/Pendulum Angle Watchdog: Y' */
  P_SRV02XAngleWatchdog_Experim_T SRV02XAngleWatchdog_n;/* '<S6>/SRV02 X Angle Watchdog' */
  P_PendulumAngleWatchdogX_Expe_T PendulumAngleWatchdogX_d;/* '<S6>/Pendulum Angle Watchdog: X' */
};

/* Real-time Model Data Structure */
struct tag_RTM_ExperimentRobustContr_T {
  const char_T *path;
  const char_T *modelName;
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;
  void *blockIO;
  const void *constBlockIO;
  void *defaultParam;
  ZCSigState *prevZCSigState;
  real_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  void *zcSignalValues;
  void *inputs;
  void *outputs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T blkStateChange;
  real_T odeF[1][64];
  ODE1_IntgData intgData;
  void *dwork;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
    uint32_T options;
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
    void *xpcData;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    void *timingData;
    real_T *varNextHitTimesList;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_ExperimentRobustControl_T ExperimentRobustControl_P;

/* Block signals (auto storage) */
extern B_ExperimentRobustControl_T ExperimentRobustControl_B;

/* Continuous states (auto storage) */
extern X_ExperimentRobustControl_T ExperimentRobustControl_X;

/* Block states (auto storage) */
extern DW_ExperimentRobustControl_T ExperimentRobustControl_DW;

/* Model entry point functions */
extern void ExperimentRobustControl_initialize(void);
extern void ExperimentRobustControl_output(void);
extern void ExperimentRobustControl_update(void);
extern void ExperimentRobustControl_terminate(void);

/*====================*
 * External functions *
 *====================*/
extern ExperimentRobustControl_rtModel *ExperimentRobustControl(void);
extern void MdlInitializeSizes(void);
extern void MdlInitializeSampleTimes(void);
extern void MdlInitialize(void);
extern void MdlStart(void);
extern void MdlOutputs(int_T tid);
extern void MdlUpdate(int_T tid);
extern void MdlTerminate(void);

/* Real-time Model object */
extern RT_MODEL_ExperimentRobustCont_T *const ExperimentRobustControl_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'ExperimentRobustControl'
 * '<S1>'   : 'ExperimentRobustControl/H Infinty controller'
 * '<S2>'   : 'ExperimentRobustControl/LQG Labb'
 * '<S3>'   : 'ExperimentRobustControl/SRV02+2DIP-E'
 * '<S4>'   : 'ExperimentRobustControl/LQG Labb/Kalman Estimator: X'
 * '<S5>'   : 'ExperimentRobustControl/LQG Labb/Kalman Estimator: Y'
 * '<S6>'   : 'ExperimentRobustControl/LQG Labb/Make State X_x'
 * '<S7>'   : 'ExperimentRobustControl/LQG Labb/Make State X_y'
 * '<S8>'   : 'ExperimentRobustControl/LQG Labb/Radians to Degrees'
 * '<S9>'   : 'ExperimentRobustControl/LQG Labb/Make State X_x/Pendulum Angle Watchdog: X'
 * '<S10>'  : 'ExperimentRobustControl/LQG Labb/Make State X_x/SRV02 X Angle Watchdog'
 * '<S11>'  : 'ExperimentRobustControl/LQG Labb/Make State X_x/Pendulum Angle Watchdog: X/Stop with Message'
 * '<S12>'  : 'ExperimentRobustControl/LQG Labb/Make State X_x/Pendulum Angle Watchdog: X/Stop with Message/Compare'
 * '<S13>'  : 'ExperimentRobustControl/LQG Labb/Make State X_x/SRV02 X Angle Watchdog/Stop with Message'
 * '<S14>'  : 'ExperimentRobustControl/LQG Labb/Make State X_x/SRV02 X Angle Watchdog/Stop with Message/Compare'
 * '<S15>'  : 'ExperimentRobustControl/LQG Labb/Make State X_y/Pendulum Angle Watchdog: Y'
 * '<S16>'  : 'ExperimentRobustControl/LQG Labb/Make State X_y/SRV02 Y Angle Watchdog'
 * '<S17>'  : 'ExperimentRobustControl/LQG Labb/Make State X_y/Pendulum Angle Watchdog: Y/Stop with Message'
 * '<S18>'  : 'ExperimentRobustControl/LQG Labb/Make State X_y/Pendulum Angle Watchdog: Y/Stop with Message/Compare'
 * '<S19>'  : 'ExperimentRobustControl/LQG Labb/Make State X_y/SRV02 Y Angle Watchdog/Stop with Message'
 * '<S20>'  : 'ExperimentRobustControl/LQG Labb/Make State X_y/SRV02 Y Angle Watchdog/Stop with Message/Compare'
 * '<S21>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Actuator Electrical  Dynamics'
 * '<S22>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Compare To Constant'
 * '<S23>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_x'
 * '<S24>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_y'
 * '<S25>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Reset Encoders'
 * '<S26>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_x/Pendulum Angle Watchdog: X'
 * '<S27>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_x/SRV02 X Angle Watchdog'
 * '<S28>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_x/Pendulum Angle Watchdog: X/Stop with Message'
 * '<S29>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_x/Pendulum Angle Watchdog: X/Stop with Message/Compare'
 * '<S30>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_x/SRV02 X Angle Watchdog/Stop with Message'
 * '<S31>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_x/SRV02 X Angle Watchdog/Stop with Message/Compare'
 * '<S32>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_y/Pendulum Angle Watchdog: Y'
 * '<S33>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_y/SRV02 Y Angle Watchdog'
 * '<S34>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_y/Pendulum Angle Watchdog: Y/Stop with Message'
 * '<S35>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_y/Pendulum Angle Watchdog: Y/Stop with Message/Compare'
 * '<S36>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_y/SRV02 Y Angle Watchdog/Stop with Message'
 * '<S37>'  : 'ExperimentRobustControl/SRV02+2DIP-E/Make State X_y/SRV02 Y Angle Watchdog/Stop with Message/Compare'
 */
#endif                                 /* RTW_HEADER_ExperimentRobustControl_h_ */
