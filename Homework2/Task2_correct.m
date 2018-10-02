clear all; clc;

%% CHANGE CALC OF H AND MAKE PATTERN PICKING STOCHASTIC

input_raw = load('input_data_numeric.csv');
x = input_raw(:,2:end);
T = 10^5;

eta = 0.02;

ABCDEF = [1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1; ...
    1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1; ...
    1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1; ...
    1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1, 1, 1, 1, -1, 1; ...
    -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1; ...
    1, 1, 1, -1, -1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1];

g = @(b)(tanh(b));
g_prim = @(b)(sech(b)^2);
answer = zeros(size(ABCDEF,1),1);
for yep = 1:size(ABCDEF,1)
    A = ABCDEF(yep,:);
    H_min = 100000;
    for attempt = 1:10
        w = rand(4,1)*0.4 - 0.2;
        theta = rand(1)*2 - 1;
        O = g(0.5*(w'*x' - theta));
        for trial = 1:T
            pattern = randi(16);
            V = x(pattern,:);
            O(pattern) = g(0.5*(w'*V' - theta));
            
            delta_out = g_prim(0.5*(w'*V' - theta))*(A(pattern)-O(pattern));
            
            %         update_weight = randi(4);
            delta_w = eta*delta_out*V';
            w = w + delta_w;
            
            delta_theta = -eta*delta_out;
            theta = theta + delta_theta;
            H = (1/2)*(A - O)*(A - O)'; %H = (1/2)*sum((A - O).^2);
            if H<0.001
                answer(yep, 1) = H;
                break
            end
            if H<H_min
                H_min = H;
            end
        end
        if H<0.001
            answer(yep, 1) = H;
            break
        end
    end
    if answer(yep,1) == 0
        answer(yep,1) = H_min;
    end
    yep
end