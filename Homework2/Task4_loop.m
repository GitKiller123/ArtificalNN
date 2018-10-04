clear all; clc;

training_set = load('training_set.csv');
training_set_in = training_set(:,1:2);
training_set_out = training_set(:,3);

validation_set = load('validation_set.csv');
validation_set_in = validation_set(:,1:2);
validation_set_out = validation_set(:,3);

T = 10^5;
eta = 0.02;
size_batch = 20;
M1 = 3; %Amount of neurons in 1st hidden layer
M2 = 2; %Amount of neurons in 2nd hidden layer

g = @(b)(tanh(b));  %Activation function
g_prim = @(b)(sech(b)^2); %Derivative of the activation function
b = @(w,V,theta,i,l)(0.5*(w{l}(:,i)'*V{l-1}' - theta{l}(i)));
delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l-1)));

N_neurons = {2, M1, M2, 1}; %Amount of Neurons in each layer

layers = length(N_neurons); %Amount of layers
% for t = 1:100
for l = 2:length(N_neurons)
    w{l} = rand(N_neurons{l-1},N_neurons{l})*0.4-0.2; %Creating starting weights for each layer
    theta{l} = rand(N_neurons{l},1)*2-1; %Creating starting thresholds for each layer
end

b = 0;
for trials = 1:10^5
    pattern = randi(size(training_set_in,1));
    V{1} = training_set_in(pattern,:);
    for l = 2:layers
        for k = 1:N_neurons{l}
            for j = 1:N_neurons{l-1}
                b = b + w{l}(j,k)*V{l-1}(j)-theta{l}(k);
            end
            V{l}(k) = tanh(b);
            b = 0;
        end
    end
    for k = 1:N_neurons{layers}
        for j = 1:N_neurons{layers-1}
            b = b + w{layers}(j,k)*V{layers-1}(j)-theta{layers}(k);
        end
        delta{layers}(k) = sech(b)^2*(training_set_out(pattern,k)-V{layers}(k));
        b = 0;
    end
    delta_temp = 0;
    for l = layers:-1:3
        for j = 1:N_neurons{l-1}
            for k = 1:N_neurons{l}
                for i = 1:N_neurons{l-2}
                    b = b + w{l-1}(i,k)*V{l-2}(i)-theta{l-1}(k);
                end
                delta_temp = delta_temp + delta{l}(k)*w{l}(j,k)*sech(b)^2;
            end
            delta{l-1}(j) = delta_temp;
            delta_temp = 0;
        end
    end
    for l = 2:layers
        for i = 1:N_neurons{l}
            for j = 1:N_neurons{l-1}
                dw(j,i) = eta*delta{l}(i)*V{l-1}(j);
            end
            dtheta(i) = -eta*delta{l}(i);
        end
        w{l} = w{l} + dw;
        theta{l} = theta{l} + dtheta;
        dw = [];
        dtheta = [];
    end
end
%% VALIDATION
for pattern = 1:size(validation_set_in,1)
    V{1} = validation_set_in(pattern,:);
    for l = 2:layers
        for k = 1:N_neurons{l}
            for j = 1:N_neurons{l-1}
                b = b + w{l}(j,k)*V{l-1}(j)-theta{l}(k);
            end
            V{l}(k) = tanh(b);
            b = 0;
        end
    end
    O(pattern,:) = V{end}(:);
end
C = (1/(2*length(validation_set_out)))*sum(abs(sign(O)-validation_set_out))
    