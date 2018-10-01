clear all; clc;

input_raw = load('input_data_numeric.csv');
x = input_raw(:,2:end);
T = 10^5;

w = rand(4,1)*0.4 - 0.2;
theta = rand(1)*2 - 1;
eta = 0.02;

A = [1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1];

g = @(b)(tanh(b));
g_prim = @(b)(sech(b)^2);
for trial = 1:T
    pattern = mod(trial,length(A));
    if pattern == 0
        pattern = length(A);
    end
    
    V = x(pattern,:);
    O(pattern) = g(0.5*(w'*V' - length(V)*theta));
    
    delta_out(size(O,1),pattern) = g_prim(0.5*(w'*V' - length(V)*theta));
    
    update_weight = randi(4);
    delta_w = eta*delta_out(pattern)*V(update_weight);
    w = w + delta_w;

    delta_theta = -eta*delta_out(pattern);
    theta = theta + delta_theta;
    if pattern == length(A)
        H = (1/2)*sum(A - O)^2; %H = (1/2)*sum((A - O).^2);
        if H<100
            H
        end
    end
end