
warning off
format compact
clc
clear all
close all

colorOrder = get(gca,'colorOrder');
%% SETUP_SRV02_EXP14_2D_GANTRY
    %
    % Sets the necessary parameters to run the SRV02 Experiment #13: Position
    % Control of 2-DOF Robot laboratory using the "s_srv02_2d_robot" and
    % "q_srv02_2d_robot" Simulink diagrams.
    % 
    % Copyright (C) 2008 Quanser Consulting Inc.
    %
    clear all;
    qc_get_step_size =1/1000; %0.001;
    deltaT = qc_get_step_size;
    Ts = qc_get_step_size;
    nPoints = (50 * 3000) -1;
%
%% Initialization Settings(DONT CHANGE ANYTHING !!)
    EXT_GEAR_CONFIG = 'HIGH';
    ENCODER_TYPE = 'E';
    TACH_OPTION = 'YES';
    LOAD_TYPE = 'NONE';
    K_AMP = 1;
    AMP_TYPE = 'VoltPAQ';
    VMAX_DAC = 10;
    ROTPEN_OPTION = '2DGANTRY-E';
    PEND_TYPE = 'MEDIUM_12IN';
    THETA_MAX = 35 * pi/180;
    ALPHA_MAX = 15.0 * pi/180;
    CONTROL_TYPE = 'AUTO';   
    X0 = pi/180*[0, 0, 0, 0];
    [ Rm, kt, km, Kg, eta_g, Beq, Jm, Jeq_noload, eta_m, K_POT, K_TACH, K_ENC, VMAX_AMP, IMAX_AMP ] = config_srv02( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, AMP_TYPE, LOAD_TYPE );
    [ g, mp, Lp, lp, Jp_cm, Bp, RtpnOp, RtpnOff, K_POT_PEN ] = config_sp( PEND_TYPE, ROTPEN_OPTION );
    [ Lb, Jarm, K_POT_2DP, K_ENC_2DP ] = config_2d_gantry( Jeq_noload );
    K_ENC_2DIP = [-1,1].*K_ENC_2DP;
    wcf_1 = 2 * pi * 5;
    zetaf_1 = 0.9;
    wcf_2 = wcf_1;
    zetaf_2 = zetaf_1;
%
%%%
%%%%
%%%%%
%%%%%% DO NOT CHANGE ANYTHING ABOVE THIS AREA !! Place your code below.
%%%%
%%%
%

%%















%%
%
%%%
%%%%
%%%%%
%%%%%% Closed loop simulation environment
%%%%
%%%
%

[ah,bh,ch,dh] = ssdata(Chinf);

figure(1)
clf;
simTime = 10;
xinit=(pi/180)*[3 3 0 0 3 3 0 0];

try
sim('Simhinf.slx')
subplot(2,1,1)
hold on
    plot(simStates.Time,simStates.Data(:,1),'linewidth',2)
    plot(simStates.Time,simStates.Data(:,3),'linewidth',2)    
    plot(simStates.Time,simStates.Data(:,2),'--','linewidth',2)
    plot(simStates.Time,simStates.Data(:,4),'--','linewidth',2)
    legend('thetaX','alphaX','thetaY','alphaY','NthetaX','NalphaX','NthetaY','NalphaY')
    ylabel('Angle [DEG]')

subplot(2,1,2)
hold on
    plot(simVoltage.Time, simVoltage.Data(:,1),'linewidth',2)
    plot(simVoltage.Time, simVoltage.Data(:,2),'--','linewidth' ,2)
    legend('Voltage X','Voltage Y')
    ylabel('Voltage [V]')
    axis([0 simTime -11 11])
catch e
    disp('Simulation failed')
end
OCL=1;

%%
%load('labb_LQR_controller.mat')
%open('ExperimentRobustControl')
