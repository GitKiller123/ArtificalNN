  function targMap = targDataMap(),

  ;%***********************
  ;% Create Parameter Map *
  ;%***********************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 23;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc paramMap
    ;%
    paramMap.nSections           = nTotSects;
    paramMap.sectIdxOffset       = sectIdxOffset;
      paramMap.sections(nTotSects) = dumSection; %prealloc
    paramMap.nTotData            = -1;
    
    ;%
    ;% Auto data (ExperimentRobustControl_P)
    ;%
      section.nData     = 31;
      section.data(31)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.ALPHA_MAX
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_P.K_AMP
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% ExperimentRobustControl_P.K_CT_labb
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% ExperimentRobustControl_P.K_ENC
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 7;
	
	  ;% ExperimentRobustControl_P.K_ENC_2DIP
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 8;
	
	  ;% ExperimentRobustControl_P.L_CT_labb
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 10;
	
	  ;% ExperimentRobustControl_P.OCL
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 18;
	
	  ;% ExperimentRobustControl_P.THETA_MAX
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 19;
	
	  ;% ExperimentRobustControl_P.VMAX_AMP
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 20;
	
	  ;% ExperimentRobustControl_P.VMAX_DAC
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 21;
	
	  ;% ExperimentRobustControl_P.ah
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 22;
	
	  ;% ExperimentRobustControl_P.bh
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 2938;
	
	  ;% ExperimentRobustControl_P.ch
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 3154;
	
	  ;% ExperimentRobustControl_P.HILInitialize_analog_input_maxi
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 3262;
	
	  ;% ExperimentRobustControl_P.HILInitialize_analog_input_mini
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 3263;
	
	  ;% ExperimentRobustControl_P.HILInitialize_analog_output_max
	  section.data(16).logicalSrcIdx = 15;
	  section.data(16).dtTransOffset = 3264;
	
	  ;% ExperimentRobustControl_P.HILInitialize_analog_output_min
	  section.data(17).logicalSrcIdx = 16;
	  section.data(17).dtTransOffset = 3265;
	
	  ;% ExperimentRobustControl_P.CompareToConstant_const
	  section.data(18).logicalSrcIdx = 17;
	  section.data(18).dtTransOffset = 3266;
	
	  ;% ExperimentRobustControl_P.HILInitialize_final_analog_outp
	  section.data(19).logicalSrcIdx = 18;
	  section.data(19).dtTransOffset = 3267;
	
	  ;% ExperimentRobustControl_P.HILInitialize_final_pwm_outputs
	  section.data(20).logicalSrcIdx = 19;
	  section.data(20).dtTransOffset = 3268;
	
	  ;% ExperimentRobustControl_P.HILInitialize_initial_analog_ou
	  section.data(21).logicalSrcIdx = 20;
	  section.data(21).dtTransOffset = 3269;
	
	  ;% ExperimentRobustControl_P.HILInitialize_initial_pwm_outpu
	  section.data(22).logicalSrcIdx = 21;
	  section.data(22).dtTransOffset = 3270;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_frequency
	  section.data(23).logicalSrcIdx = 22;
	  section.data(23).dtTransOffset = 3271;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_leading_deadb
	  section.data(24).logicalSrcIdx = 23;
	  section.data(24).dtTransOffset = 3272;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_trailing_dead
	  section.data(25).logicalSrcIdx = 24;
	  section.data(25).dtTransOffset = 3273;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_other_outputs
	  section.data(26).logicalSrcIdx = 25;
	  section.data(26).dtTransOffset = 3274;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_other_outpu_c
	  section.data(27).logicalSrcIdx = 26;
	  section.data(27).dtTransOffset = 3275;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_other_outpu_p
	  section.data(28).logicalSrcIdx = 27;
	  section.data(28).dtTransOffset = 3276;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_other_outpu_h
	  section.data(29).logicalSrcIdx = 28;
	  section.data(29).dtTransOffset = 3277;
	
	  ;% ExperimentRobustControl_P.HILInitialize_watchdog_analog_o
	  section.data(30).logicalSrcIdx = 29;
	  section.data(30).dtTransOffset = 3278;
	
	  ;% ExperimentRobustControl_P.HILInitialize_watchdog_pwm_outp
	  section.data(31).logicalSrcIdx = 30;
	  section.data(31).dtTransOffset = 3279;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
      section.nData     = 9;
      section.data(9)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.HILReadEncoderTimebase_clock
	  section.data(1).logicalSrcIdx = 31;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_P.HILSetEncoderCounts_counts
	  section.data(2).logicalSrcIdx = 32;
	  section.data(2).dtTransOffset = 1;
	
	  ;% ExperimentRobustControl_P.HILInitialize_hardware_clocks
	  section.data(3).logicalSrcIdx = 33;
	  section.data(3).dtTransOffset = 2;
	
	  ;% ExperimentRobustControl_P.HILInitialize_initial_encoder_c
	  section.data(4).logicalSrcIdx = 34;
	  section.data(4).dtTransOffset = 5;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_alignment
	  section.data(5).logicalSrcIdx = 35;
	  section.data(5).dtTransOffset = 6;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_configuration
	  section.data(6).logicalSrcIdx = 36;
	  section.data(6).dtTransOffset = 7;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_modes
	  section.data(7).logicalSrcIdx = 37;
	  section.data(7).dtTransOffset = 8;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_polarity
	  section.data(8).logicalSrcIdx = 38;
	  section.data(8).dtTransOffset = 9;
	
	  ;% ExperimentRobustControl_P.HILInitialize_watchdog_digital_
	  section.data(9).logicalSrcIdx = 39;
	  section.data(9).dtTransOffset = 10;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(2) = section;
      clear section
      
      section.nData     = 11;
      section.data(11)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.HILInitialize_analog_input_chan
	  section.data(1).logicalSrcIdx = 40;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_P.HILInitialize_analog_output_cha
	  section.data(2).logicalSrcIdx = 41;
	  section.data(2).dtTransOffset = 8;
	
	  ;% ExperimentRobustControl_P.HILSetEncoderCounts_channels
	  section.data(3).logicalSrcIdx = 42;
	  section.data(3).dtTransOffset = 16;
	
	  ;% ExperimentRobustControl_P.HILReadEncoderTimebase_channels
	  section.data(4).logicalSrcIdx = 43;
	  section.data(4).dtTransOffset = 20;
	
	  ;% ExperimentRobustControl_P.HILWriteAnalog_channels
	  section.data(5).logicalSrcIdx = 44;
	  section.data(5).dtTransOffset = 24;
	
	  ;% ExperimentRobustControl_P.HILWriteDigital_channels
	  section.data(6).logicalSrcIdx = 45;
	  section.data(6).dtTransOffset = 26;
	
	  ;% ExperimentRobustControl_P.HILInitialize_digital_output_ch
	  section.data(7).logicalSrcIdx = 46;
	  section.data(7).dtTransOffset = 30;
	
	  ;% ExperimentRobustControl_P.HILInitialize_encoder_channels
	  section.data(8).logicalSrcIdx = 47;
	  section.data(8).dtTransOffset = 38;
	
	  ;% ExperimentRobustControl_P.HILInitialize_pwm_channels
	  section.data(9).logicalSrcIdx = 48;
	  section.data(9).dtTransOffset = 46;
	
	  ;% ExperimentRobustControl_P.HILInitialize_quadrature
	  section.data(10).logicalSrcIdx = 49;
	  section.data(10).dtTransOffset = 54;
	
	  ;% ExperimentRobustControl_P.HILReadEncoderTimebase_samples_
	  section.data(11).logicalSrcIdx = 50;
	  section.data(11).dtTransOffset = 55;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(3) = section;
      clear section
      
      section.nData     = 35;
      section.data(35)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.HILInitialize_active
	  section.data(1).logicalSrcIdx = 51;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_P.HILInitialize_final_digital_out
	  section.data(2).logicalSrcIdx = 52;
	  section.data(2).dtTransOffset = 1;
	
	  ;% ExperimentRobustControl_P.HILInitialize_initial_digital_o
	  section.data(3).logicalSrcIdx = 53;
	  section.data(3).dtTransOffset = 2;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_input_
	  section.data(4).logicalSrcIdx = 54;
	  section.data(4).dtTransOffset = 3;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_inpu_o
	  section.data(5).logicalSrcIdx = 55;
	  section.data(5).dtTransOffset = 4;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_output
	  section.data(6).logicalSrcIdx = 56;
	  section.data(6).dtTransOffset = 5;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_outp_d
	  section.data(7).logicalSrcIdx = 57;
	  section.data(7).dtTransOffset = 6;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_outp_j
	  section.data(8).logicalSrcIdx = 58;
	  section.data(8).dtTransOffset = 7;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_outp_m
	  section.data(9).logicalSrcIdx = 59;
	  section.data(9).dtTransOffset = 8;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_out_mi
	  section.data(10).logicalSrcIdx = 60;
	  section.data(10).dtTransOffset = 9;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_outp_n
	  section.data(11).logicalSrcIdx = 61;
	  section.data(11).dtTransOffset = 10;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_analog_outp_k
	  section.data(12).logicalSrcIdx = 62;
	  section.data(12).dtTransOffset = 11;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_clock_frequen
	  section.data(13).logicalSrcIdx = 63;
	  section.data(13).dtTransOffset = 12;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_clock_frequ_l
	  section.data(14).logicalSrcIdx = 64;
	  section.data(14).dtTransOffset = 13;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_clock_params_
	  section.data(15).logicalSrcIdx = 65;
	  section.data(15).dtTransOffset = 14;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_clock_param_o
	  section.data(16).logicalSrcIdx = 66;
	  section.data(16).dtTransOffset = 15;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_digital_outpu
	  section.data(17).logicalSrcIdx = 67;
	  section.data(17).dtTransOffset = 16;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_digital_out_j
	  section.data(18).logicalSrcIdx = 68;
	  section.data(18).dtTransOffset = 17;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_digital_out_a
	  section.data(19).logicalSrcIdx = 69;
	  section.data(19).dtTransOffset = 18;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_digital_out_f
	  section.data(20).logicalSrcIdx = 70;
	  section.data(20).dtTransOffset = 19;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_digital_out_b
	  section.data(21).logicalSrcIdx = 71;
	  section.data(21).dtTransOffset = 20;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_digital_ou_bw
	  section.data(22).logicalSrcIdx = 72;
	  section.data(22).dtTransOffset = 21;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_digital_ou_j0
	  section.data(23).logicalSrcIdx = 73;
	  section.data(23).dtTransOffset = 22;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_encoder_count
	  section.data(24).logicalSrcIdx = 74;
	  section.data(24).dtTransOffset = 23;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_encoder_cou_e
	  section.data(25).logicalSrcIdx = 75;
	  section.data(25).dtTransOffset = 24;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_encoder_param
	  section.data(26).logicalSrcIdx = 76;
	  section.data(26).dtTransOffset = 25;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_encoder_par_m
	  section.data(27).logicalSrcIdx = 77;
	  section.data(27).dtTransOffset = 26;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_other_outpu_i
	  section.data(28).logicalSrcIdx = 78;
	  section.data(28).dtTransOffset = 27;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_a
	  section.data(29).logicalSrcIdx = 79;
	  section.data(29).dtTransOffset = 28;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_l
	  section.data(30).logicalSrcIdx = 80;
	  section.data(30).dtTransOffset = 29;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_n
	  section.data(31).logicalSrcIdx = 81;
	  section.data(31).dtTransOffset = 30;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_c
	  section.data(32).logicalSrcIdx = 82;
	  section.data(32).dtTransOffset = 31;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_pwm_outputs_o
	  section.data(33).logicalSrcIdx = 83;
	  section.data(33).dtTransOffset = 32;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_pwm_params_at
	  section.data(34).logicalSrcIdx = 84;
	  section.data(34).dtTransOffset = 33;
	
	  ;% ExperimentRobustControl_P.HILInitialize_set_pwm_params__m
	  section.data(35).logicalSrcIdx = 85;
	  section.data(35).dtTransOffset = 34;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(4) = section;
      clear section
      
      section.nData     = 16;
      section.data(16)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.StateSpace_X0
	  section.data(1).logicalSrcIdx = 86;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_P.Integrator_IC
	  section.data(2).logicalSrcIdx = 87;
	  section.data(2).dtTransOffset = 1;
	
	  ;% ExperimentRobustControl_P.Integrator_IC_b
	  section.data(3).logicalSrcIdx = 88;
	  section.data(3).dtTransOffset = 2;
	
	  ;% ExperimentRobustControl_P.Integrator_IC_p
	  section.data(4).logicalSrcIdx = 89;
	  section.data(4).dtTransOffset = 3;
	
	  ;% ExperimentRobustControl_P.Integrator_IC_px
	  section.data(5).logicalSrcIdx = 90;
	  section.data(5).dtTransOffset = 4;
	
	  ;% ExperimentRobustControl_P.DirectionConventionRightHandsys
	  section.data(6).logicalSrcIdx = 91;
	  section.data(6).dtTransOffset = 5;
	
	  ;% ExperimentRobustControl_P.EnableVoltPAQX2X4_Value
	  section.data(7).logicalSrcIdx = 92;
	  section.data(7).dtTransOffset = 6;
	
	  ;% ExperimentRobustControl_P.Gain2_Gain
	  section.data(8).logicalSrcIdx = 93;
	  section.data(8).dtTransOffset = 10;
	
	  ;% ExperimentRobustControl_P.Gain_Gain
	  section.data(9).logicalSrcIdx = 94;
	  section.data(9).dtTransOffset = 11;
	
	  ;% ExperimentRobustControl_P.Gain_Gain_d
	  section.data(10).logicalSrcIdx = 95;
	  section.data(10).dtTransOffset = 12;
	
	  ;% ExperimentRobustControl_P.syscc_x_Gain
	  section.data(11).logicalSrcIdx = 96;
	  section.data(11).dtTransOffset = 13;
	
	  ;% ExperimentRobustControl_P.syscb_x_Gain
	  section.data(12).logicalSrcIdx = 97;
	  section.data(12).dtTransOffset = 21;
	
	  ;% ExperimentRobustControl_P.sysca_x_Gain
	  section.data(13).logicalSrcIdx = 98;
	  section.data(13).dtTransOffset = 25;
	
	  ;% ExperimentRobustControl_P.sysdc_x_Gain
	  section.data(14).logicalSrcIdx = 99;
	  section.data(14).dtTransOffset = 41;
	
	  ;% ExperimentRobustControl_P.sysdb_x_Gain
	  section.data(15).logicalSrcIdx = 100;
	  section.data(15).dtTransOffset = 49;
	
	  ;% ExperimentRobustControl_P.sysda_x_Gain
	  section.data(16).logicalSrcIdx = 101;
	  section.data(16).dtTransOffset = 53;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(5) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.ManualSwitch2_CurrentSetting
	  section.data(1).logicalSrcIdx = 102;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(6) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.HILSetEncoderCounts_Active
	  section.data(1).logicalSrcIdx = 103;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_P.HILReadEncoderTimebase_Active
	  section.data(2).logicalSrcIdx = 104;
	  section.data(2).dtTransOffset = 1;
	
	  ;% ExperimentRobustControl_P.HILWriteAnalog_Active
	  section.data(3).logicalSrcIdx = 105;
	  section.data(3).dtTransOffset = 2;
	
	  ;% ExperimentRobustControl_P.HILWriteDigital_Active
	  section.data(4).logicalSrcIdx = 106;
	  section.data(4).dtTransOffset = 3;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(7) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02YAngleWatchdog.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 107;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(8) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02YAngleWatchdog.Constant_Value
	  section.data(1).logicalSrcIdx = 108;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(9) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogY.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 109;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(10) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogY.Constant_Value
	  section.data(1).logicalSrcIdx = 110;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(11) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02XAngleWatchdog.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 111;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(12) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02XAngleWatchdog.Constant_Value
	  section.data(1).logicalSrcIdx = 112;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(13) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogX.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 113;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(14) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogX.Constant_Value
	  section.data(1).logicalSrcIdx = 114;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(15) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02YAngleWatchdog_l.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 115;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(16) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02YAngleWatchdog_l.Constant_Value
	  section.data(1).logicalSrcIdx = 116;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(17) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogY_a.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 117;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(18) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogY_a.Constant_Value
	  section.data(1).logicalSrcIdx = 118;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(19) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02XAngleWatchdog_n.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 119;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(20) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.SRV02XAngleWatchdog_n.Constant_Value
	  section.data(1).logicalSrcIdx = 120;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(21) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogX_d.StopwithMessage_message_icon
	  section.data(1).logicalSrcIdx = 121;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(22) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_P.PendulumAngleWatchdogX_d.Constant_Value
	  section.data(1).logicalSrcIdx = 122;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(23) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (parameter)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    paramMap.nTotData = nTotData;
    


  ;%**************************
  ;% Create Block Output Map *
  ;%**************************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 18;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc sigMap
    ;%
    sigMap.nSections           = nTotSects;
    sigMap.sectIdxOffset       = sectIdxOffset;
      sigMap.sections(nTotSects) = dumSection; %prealloc
    sigMap.nTotData            = -1;
    
    ;%
    ;% Auto data (ExperimentRobustControl_B)
    ;%
      section.nData     = 41;
      section.data(41)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.HILReadEncoderTimebase_o1
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.HILReadEncoderTimebase_o2
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% ExperimentRobustControl_B.HILReadEncoderTimebase_o3
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% ExperimentRobustControl_B.HILReadEncoderTimebase_o4
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% ExperimentRobustControl_B.uDIPEncoderCalibrationradcount
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% ExperimentRobustControl_B.EncoderCalibrationradcount
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% ExperimentRobustControl_B.uDIPEncoderYCalibrationradcount
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% ExperimentRobustControl_B.EncoderCalibrationYradcount
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 7;
	
	  ;% ExperimentRobustControl_B.StateSpace
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 8;
	
	  ;% ExperimentRobustControl_B.Integrator
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 10;
	
	  ;% ExperimentRobustControl_B.Integrator_m
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 14;
	
	  ;% ExperimentRobustControl_B.Integrator_k
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 15;
	
	  ;% ExperimentRobustControl_B.Integrator_c
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 19;
	
	  ;% ExperimentRobustControl_B.TmpSignalConversionAtLQRGainYIn
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 20;
	
	  ;% ExperimentRobustControl_B.ManualSwitch2
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 25;
	
	  ;% ExperimentRobustControl_B.Gain1
	  section.data(16).logicalSrcIdx = 15;
	  section.data(16).dtTransOffset = 27;
	
	  ;% ExperimentRobustControl_B.DirectionConventionRightHandsys
	  section.data(17).logicalSrcIdx = 16;
	  section.data(17).dtTransOffset = 29;
	
	  ;% ExperimentRobustControl_B.AmplifierSaturationV
	  section.data(18).logicalSrcIdx = 17;
	  section.data(18).dtTransOffset = 31;
	
	  ;% ExperimentRobustControl_B.InverseAmplifierGainVV
	  section.data(19).logicalSrcIdx = 18;
	  section.data(19).dtTransOffset = 33;
	
	  ;% ExperimentRobustControl_B.DACBSaturationV
	  section.data(20).logicalSrcIdx = 19;
	  section.data(20).dtTransOffset = 35;
	
	  ;% ExperimentRobustControl_B.AmplifierGainVV
	  section.data(21).logicalSrcIdx = 20;
	  section.data(21).dtTransOffset = 37;
	
	  ;% ExperimentRobustControl_B.Clock
	  section.data(22).logicalSrcIdx = 21;
	  section.data(22).dtTransOffset = 39;
	
	  ;% ExperimentRobustControl_B.Gain2
	  section.data(23).logicalSrcIdx = 22;
	  section.data(23).dtTransOffset = 40;
	
	  ;% ExperimentRobustControl_B.Gain
	  section.data(24).logicalSrcIdx = 23;
	  section.data(24).dtTransOffset = 42;
	
	  ;% ExperimentRobustControl_B.TmpSignalConversionAtStateSpace
	  section.data(25).logicalSrcIdx = 24;
	  section.data(25).dtTransOffset = 44;
	
	  ;% ExperimentRobustControl_B.Gain_l
	  section.data(26).logicalSrcIdx = 25;
	  section.data(26).dtTransOffset = 48;
	
	  ;% ExperimentRobustControl_B.syscc_x
	  section.data(27).logicalSrcIdx = 26;
	  section.data(27).dtTransOffset = 53;
	
	  ;% ExperimentRobustControl_B.Sum1
	  section.data(28).logicalSrcIdx = 27;
	  section.data(28).dtTransOffset = 55;
	
	  ;% ExperimentRobustControl_B.KalmanFilterGain
	  section.data(29).logicalSrcIdx = 28;
	  section.data(29).dtTransOffset = 57;
	
	  ;% ExperimentRobustControl_B.syscb_x
	  section.data(30).logicalSrcIdx = 29;
	  section.data(30).dtTransOffset = 61;
	
	  ;% ExperimentRobustControl_B.sysca_x
	  section.data(31).logicalSrcIdx = 30;
	  section.data(31).dtTransOffset = 65;
	
	  ;% ExperimentRobustControl_B.Sum
	  section.data(32).logicalSrcIdx = 31;
	  section.data(32).dtTransOffset = 69;
	
	  ;% ExperimentRobustControl_B.sysdc_x
	  section.data(33).logicalSrcIdx = 32;
	  section.data(33).dtTransOffset = 73;
	
	  ;% ExperimentRobustControl_B.Sum1_k
	  section.data(34).logicalSrcIdx = 33;
	  section.data(34).dtTransOffset = 75;
	
	  ;% ExperimentRobustControl_B.KalmanFilterGain_e
	  section.data(35).logicalSrcIdx = 34;
	  section.data(35).dtTransOffset = 77;
	
	  ;% ExperimentRobustControl_B.sysdb_x
	  section.data(36).logicalSrcIdx = 35;
	  section.data(36).dtTransOffset = 81;
	
	  ;% ExperimentRobustControl_B.sysda_x
	  section.data(37).logicalSrcIdx = 36;
	  section.data(37).dtTransOffset = 85;
	
	  ;% ExperimentRobustControl_B.Sum_j
	  section.data(38).logicalSrcIdx = 37;
	  section.data(38).dtTransOffset = 89;
	
	  ;% ExperimentRobustControl_B.TmpSignalConversionAtLQRGainXIn
	  section.data(39).logicalSrcIdx = 38;
	  section.data(39).dtTransOffset = 93;
	
	  ;% ExperimentRobustControl_B.LQRGainX
	  section.data(40).logicalSrcIdx = 39;
	  section.data(40).dtTransOffset = 98;
	
	  ;% ExperimentRobustControl_B.LQRGainY
	  section.data(41).logicalSrcIdx = 40;
	  section.data(41).dtTransOffset = 99;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.Compare
	  section.data(1).logicalSrcIdx = 41;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(2) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02YAngleWatchdog.Abs
	  section.data(1).logicalSrcIdx = 42;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.SRV02YAngleWatchdog.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 43;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(3) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02YAngleWatchdog.Compare
	  section.data(1).logicalSrcIdx = 44;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(4) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogY.Abs
	  section.data(1).logicalSrcIdx = 45;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogY.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 46;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(5) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogY.Compare
	  section.data(1).logicalSrcIdx = 47;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(6) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02XAngleWatchdog.Abs
	  section.data(1).logicalSrcIdx = 48;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.SRV02XAngleWatchdog.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 49;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(7) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02XAngleWatchdog.Compare
	  section.data(1).logicalSrcIdx = 50;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(8) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogX.Abs
	  section.data(1).logicalSrcIdx = 51;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogX.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 52;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(9) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogX.Compare
	  section.data(1).logicalSrcIdx = 53;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(10) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02YAngleWatchdog_l.Abs
	  section.data(1).logicalSrcIdx = 54;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.SRV02YAngleWatchdog_l.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 55;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(11) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02YAngleWatchdog_l.Compare
	  section.data(1).logicalSrcIdx = 56;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(12) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogY_a.Abs
	  section.data(1).logicalSrcIdx = 57;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogY_a.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 58;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(13) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogY_a.Compare
	  section.data(1).logicalSrcIdx = 59;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(14) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02XAngleWatchdog_n.Abs
	  section.data(1).logicalSrcIdx = 60;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.SRV02XAngleWatchdog_n.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 61;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(15) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.SRV02XAngleWatchdog_n.Compare
	  section.data(1).logicalSrcIdx = 62;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(16) = section;
      clear section
      
      section.nData     = 2;
      section.data(2)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogX_d.Abs
	  section.data(1).logicalSrcIdx = 63;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogX_d.CheckMaximumAllowedAngle
	  section.data(2).logicalSrcIdx = 64;
	  section.data(2).dtTransOffset = 1;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(17) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_B.PendulumAngleWatchdogX_d.Compare
	  section.data(1).logicalSrcIdx = 65;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(18) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (signal)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    sigMap.nTotData = nTotData;
    


  ;%*******************
  ;% Create DWork Map *
  ;%*******************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 10;
    sectIdxOffset = 18;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc dworkMap
    ;%
    dworkMap.nSections           = nTotSects;
    dworkMap.sectIdxOffset       = sectIdxOffset;
      dworkMap.sections(nTotSects) = dumSection; %prealloc
    dworkMap.nTotData            = -1;
    
    ;%
    ;% Auto data (ExperimentRobustControl_DW)
    ;%
      section.nData     = 8;
      section.data(8)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILInitialize_AIMinimums
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_AIMaximums
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 8;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_AOMinimums
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 16;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_AOMaximums
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 24;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_AOVoltages
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 32;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_FilterFrequency
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 40;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_POSortedFreqs
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 48;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_POValues
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 56;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILInitialize_Card
	  section.data(1).logicalSrcIdx = 8;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILReadEncoderTimebase_Task
	  section.data(1).logicalSrcIdx = 9;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(3) = section;
      clear section
      
      section.nData     = 5;
      section.data(5)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILWriteAnalog_PWORK
	  section.data(1).logicalSrcIdx = 10;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_DW.HILWriteDigital_PWORK
	  section.data(2).logicalSrcIdx = 11;
	  section.data(2).dtTransOffset = 1;
	
	  ;% ExperimentRobustControl_DW.Scope1_PWORK.LoggedData
	  section.data(3).logicalSrcIdx = 12;
	  section.data(3).dtTransOffset = 2;
	
	  ;% ExperimentRobustControl_DW.Scope_PWORK.LoggedData
	  section.data(4).logicalSrcIdx = 13;
	  section.data(4).dtTransOffset = 5;
	
	  ;% ExperimentRobustControl_DW.HILSetEncoderCounts_PWORK
	  section.data(5).logicalSrcIdx = 14;
	  section.data(5).dtTransOffset = 6;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(4) = section;
      clear section
      
      section.nData     = 9;
      section.data(9)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILInitialize_ClockModes
	  section.data(1).logicalSrcIdx = 15;
	  section.data(1).dtTransOffset = 0;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_DOStates
	  section.data(2).logicalSrcIdx = 16;
	  section.data(2).dtTransOffset = 3;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_QuadratureModes
	  section.data(3).logicalSrcIdx = 17;
	  section.data(3).dtTransOffset = 11;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_InitialEICounts
	  section.data(4).logicalSrcIdx = 18;
	  section.data(4).dtTransOffset = 19;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_POModeValues
	  section.data(5).logicalSrcIdx = 19;
	  section.data(5).dtTransOffset = 27;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_POAlignValues
	  section.data(6).logicalSrcIdx = 20;
	  section.data(6).dtTransOffset = 35;
	
	  ;% ExperimentRobustControl_DW.HILInitialize_POPolarityVals
	  section.data(7).logicalSrcIdx = 21;
	  section.data(7).dtTransOffset = 43;
	
	  ;% ExperimentRobustControl_DW.HILReadEncoderTimebase_Buffer
	  section.data(8).logicalSrcIdx = 22;
	  section.data(8).dtTransOffset = 51;
	
	  ;% ExperimentRobustControl_DW.HILSetEncoderCounts_Buffer
	  section.data(9).logicalSrcIdx = 23;
	  section.data(9).dtTransOffset = 55;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(5) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILInitialize_POSortedChans
	  section.data(1).logicalSrcIdx = 24;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(6) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.ResetEncoders_SubsysRanBC
	  section.data(1).logicalSrcIdx = 25;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(7) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILInitialize_DOBits
	  section.data(1).logicalSrcIdx = 26;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(8) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.HILWriteDigital_Buffer
	  section.data(1).logicalSrcIdx = 27;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(9) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% ExperimentRobustControl_DW.ResetEncoders_MODE
	  section.data(1).logicalSrcIdx = 28;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(10) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (dwork)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    dworkMap.nTotData = nTotData;
    


  ;%
  ;% Add individual maps to base struct.
  ;%

  targMap.paramMap  = paramMap;    
  targMap.signalMap = sigMap;
  targMap.dworkMap  = dworkMap;
  
  ;%
  ;% Add checksums to base struct.
  ;%


  targMap.checksum0 = 3116169259;
  targMap.checksum1 = 2044929189;
  targMap.checksum2 = 3951391620;
  targMap.checksum3 = 600318527;

