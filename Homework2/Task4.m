clear all; clc;

training_set = load('training_set.csv');
training_set_in = training_set(:,1:2);
training_set_out = training_set(:,3);
T = 10^5;
eta = 0.02;
M1 = 3; %Amount of neurons in 1st hidden layer
M2 = 2; %Amount of neurons in 2nd hidden layer

g = @(b)(tanh(b));  %Activation function
g_prim = @(b)(sech(b)^2); %Derivative of the activation function
b = @(w,V,theta,i,l)(0.5*(w{l-1}(:,i)'*V{l-1}' - theta{l-1}(i)));
delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l)));

N_neurons = {2, M1, M2, 1}; %Amount of Neurons in each layer

layers = length(N_neurons); %Amount of layers
for l = 2:length(N_neurons)
    w{l-1} = rand(N_neurons{l-1},N_neurons{l})*0.4-0.2; %Creating starting weights for each layer
    theta{l-1} = rand(N_neurons{l},1)*2-1; %Creating starting thresholds for each layer
end
%%
% TRAINING
for trial = 1:T %size(training_set_in,1)
    pattern_start = randi(size(training_set,1)-30);
    for pattern = pattern_start:pattern_start+30
        V{1} = training_set_in(pattern,:);
        for l = 2:layers
            for i = 1:N_neurons{l}
                V{l}(i) = g(b(w,V,theta,i,l));
            end
        end

        for i = 1:N_neurons{end}
            delta{layers-1}(i) = g_prim(b(w,V,theta,i,l))*(training_set_out(pattern,i)-V{end}(i));
        end
        for l = layers-1:-1:2
            for j = 1:N_neurons{l}
                for i = 1:N_neurons{l+1}
                    delta_temp(i) = delta_func(w,V,theta,delta,i,j,l);
                    %                 delta{l-1}(j) = delta{l}(i)*w{l}(j)*g_prim(b(w,V,theta,j,l));
                end
                delta{l-1}(j) = sum(delta_temp);
                delta_temp = [];
            end
        end
        for l = 1:layers-1
            delta_w = eta*delta{l}'*V{l};
%             delta_w = eta*delta{l}'*V{l};
            w{l} = w{l} + delta_w';
            
            delta_theta = -eta*delta{l}';
            theta{l} = theta{l} + delta_theta;
        end
        O_training(pattern-pattern_start+1) = V{end}(1);
    end
    training_out_curr = training_set_out(pattern_start:pattern_start+30);
    H = 0.5*(training_out_curr'-O_training)*(training_out_curr'-O_training)'
end
sum(sign(O_training')-training_set_out)

%% VALIDATION
validation_set = load('validation_set.csv');
validation_set_in = validation_set(:,1:2);
validation_set_out = validation_set(:,3);

for pattern = 1:size(validation_set_in,1)
    V{1} = validation_set_in(pattern,:);
    for l = 2:layers
        for i = 1:N_neurons{l}
            V{l}(i) = g(b(w,V,theta,i,l));
        end
    end
    O(pattern,:) = V{end}(:);
end
C = (1/(2*length(validation_set_out)))*sum(abs(sign(O)-validation_set_out));
