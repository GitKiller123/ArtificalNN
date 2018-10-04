%
%% CONFIG_2D_GANTRY
%
% Sets up and returns the 2-DOF Gantry model variables
%
% 2-DOF Gantry Module system nomenclature:
% ************************************************************************
% Input parameters:
%
% Jeq_m     Initial Equivalent Inertia of shaft without the 2-DOF robot 
%           and gantry modules                                  (kg.m^2)
%
% ************************************************************************
% Output parameters:
% Lb            Bar Linkage Length                               (m)
% Jeq         Equivalent Inertia as seen at each SRV02 Load shaft
%               with the 4-bar linkage + 2-DOF gantry assembly. (kg.m^2)
% K_POT_2DP     Potentiometer Sensitivity                       (rad/V)
% K_ENC_2DP     Encoder Resolution (in quadrature mode)         (rad/count)
%
% Copyright (C) 2008 Quanser Consulting Inc.
% Quanser Consulting Inc.
%%
function [ Lb, Jeq, K_POT_2DP, K_ENC_2DP ] = config_2d_gantry( Jeq_m )
% Potentiometer Sensitivity (rad/V)
K_POT_2DP = 14.6 * pi / 180*[-1,1];
% Encoder Resolution, for a quadrature encoder, (rad/count)
K_ENC_2DP = 2 * pi / ( 4 * 1024 )*[-1,1];
%
% Load 2D ROBOT parameters
[ Lb_in, Jeq_2drobot, Jeq_linkage ] = config_2d_robot( Jeq_m );
% Mass of 2-DOF hinge with 2 potentiometers (kg)
Mh = 0.3;
% Length of link (m)
Lb = Lb_in * 0.0254;
% Moment of inertia of linkage assembly (kg.m^2)
Jeq_assembly = Jeq_linkage + Mh * Lb^2; % = 0.0065
% Equivalent inertia with the 4-bar linkage + 2-DOF gantry assembly (kg.m^2)
Jeq = Jeq_m + Jeq_assembly; % = 0.0086