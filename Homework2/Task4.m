clear all; clc;

training_set = load('training_set.csv');
training_set_in = training_set(:,1:2);
training_set_out = training_set(:,3);
validation_set = load('validation_set.csv');
T = 10^5;
eta = 0.02;

g = @(b)(tanh(b));
g_prim = @(b)(sech(b)^2);

M1 = 3;
M2 = 2;
w1 = rand(2,M1)*0.4 - 0.2;
w2 = rand(M1,M2)*0.4 - 0.2;
theta1 = rand(M1,1)*2 - 1;
theta2 = rand(M2,1)*2 - 1;
theta_out = rand(1)*2 - 1;
    for attempt = 1:10
        w = rand(4,1)*0.4 - 0.2;
        theta = rand(1)*2 - 1;
        for trial = 1:T
            pattern = randi(size(training_set,1));
            V = training_set(pattern,:); %HIT HAR VI KÖRT
            O(pattern) = g(0.5*(w'*V' - length(V)*theta));
            
            delta_out(pattern) = g_prim(0.5*(w'*V' - length(V)*theta))*(A(pattern)-O(pattern));
            
            %         update_weight = randi(4);
            delta_w = eta*delta_out(pattern)*V';
            w = w + delta_w;
            
            delta_theta = -eta*delta_out(pattern);
            theta = theta + delta_theta;
            H = (1/2)*(A - O)*(A - O)'; %H = (1/2)*sum((A - O).^2);
            if H<0.01
                answer(yep, 1) = H;
                break
            end
            if H<H_min
                H_min = H;
            end
        end
        if H<0.01
            answer(yep, 1) = H;
            break
        end
    end