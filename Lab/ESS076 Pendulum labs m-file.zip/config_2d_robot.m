%
%% CONFIG_2DROBOT
%
% Sets up and returns the 2-DOF Manipulator model variables.
%
% 2-DOF Robot Module system nomenclature:

% ************************************************************************
% Input parameters:
% Jeq_m         Initial Equivalent Inertia of shaft without the 2-DOF   
%               robot module.                                       (kg.m^2)
%
% ************************************************************************
% Output parameters:
% Lb_in         Link Length                                         (in)
% Jeq           Equivalent Inertia as seen at each SRV02 load shaft 
%               with the 4-bar linkage.                             (kg.m^2)
% Jeq_linkage   Equivalent moment of inertia of the 4-bar linkage.  (kg.m^2)
%
% Copyright (C) 2008 Quanser Consulting Inc.
% Quanser Consulting Inc.
%
%%
function [ Lb_in, Jeq, Jeq_linkage ] = config_2d_robot( Jeq_m )
% 
% load support bar length (inch)
Lb_in = 5;
% load support bar length (m)
Lb = 5 * 0.0254;  % = 0.127
% linkage bar mass (kg)
Mb = 0.065;
% single bar inertia about its cog (kg.m^2)
Jb_cog = Mb * Lb^2 / 12; 
% single bar inertia about pivot (kg.m^2)
% Lb: distance from the pivot axis to bar cog (m)
Jb_piv = Jb_cog + Mb * (Lb/2)^2; % = 0.0011 Mb * Lb^2 / 3; 
% Jeq_linkage: equivalent inertia of the 4-bar linkage (kg.m^2)
Jeq_linkage = Jb_piv + Mb*Lb^2; % = 0.0017 old: 1.5*Jb_piv
% Jeq: equivalent inertia with the 4-bar linkage (kg.m^2)
Jeq = Jeq_m + Jeq_linkage; % = 0.0038