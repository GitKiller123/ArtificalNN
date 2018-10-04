% This function gets the disturbance vector L as functions of the nominal
% model. The parameter values are changed according to chosen pendulum_type

% Gets the disturbance vector dependent on the parameter chosen, Lp seems
% to affect the plant the most, which is what we ended up using.

function [disturbanceVector] = getDisturbanceVector(uncertaintyType,invertedPendulumParameters)

    Mp = invertedPendulumParameters.Mp;
    Lp = invertedPendulumParameters.Lp;
    Jp = invertedPendulumParameters.Jp;
    Lr = invertedPendulumParameters.Lr;
    Jr = invertedPendulumParameters.Jr;
    Dr = invertedPendulumParameters.Dr;
    Co = invertedPendulumParameters.Co;
    g = invertedPendulumParameters.g;

    switch uncertaintyType
        case 'Mp'
            alphaD = 40*Lr^2*Lp*Mp + 40*Lr^2*Lp*g - 2*Lr*Lp*Dr; % These are most likely correct. Maybe.
            thetaD = -20*Lp*Lr*g*Mp + 20*Lp*Lr*g - Lp^2 *Dr;
        case 'Lp'
            thetaD = (20*Lr*Mp^2*g*Lp* + 20*Lr*Mp^2*g) - (2*Dr*Mp*Lp + 2*Dr*Mp);
            alphaD = (20*g*Lr^2*Mp^2 + 20*Jr*g*Mp) -2*Dr*Lr*Mp;
        case 'Jp'
            thetaD = -4*Dr;
            alphaD = 0;
     end
    disturbanceVector = [ 0 0 thetaD alphaD 0]';
end