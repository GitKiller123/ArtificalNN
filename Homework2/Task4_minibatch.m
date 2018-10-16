clear all; clc;

training_set = load('training_set.csv');
training_set_in = training_set(:,1:2);
training_set_out = training_set(:,3);

validation_set = load('validation_set.csv');
validation_set_in = validation_set(:,1:2);
validation_set_out = validation_set(:,3);

T = 10^5;
mu = 0.03;
epochs = 100;
size_batch = 1;
M1 = 8; %Amount of neurons in 1st hidden layer
M2 = 5; %Amount of neurons in 2nd hidden layer

g = @(b)(tanh(b));  %Activation function
g_prim = @(b)(sech(b).^2); %Derivative of the activation function
% b = @(w,V,theta,i,l)(w{l}(:,i)'*V{l-1} - theta{l}(i));
delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l-1)));

N_neurons = {2, M1, M2, 1}; %Amount of Neurons in each layer

layers = length(N_neurons); %Amount of layers

C_best = 0.17;

for l = 2:length(N_neurons)
    w{l} = randn(N_neurons{l-1},N_neurons{l}); %Creating starting weights for each layer
    theta{l} = zeros(N_neurons{l},1); %Creating starting thresholds for each layer
end
for epochs = 1:T
    for trial = 1:500
        % TRAINING
        
        if trial + size_batch > size(training_set_in,1)
            break
        end
        for l = 1:layers
            V_tot{l} = zeros(N_neurons{l},1);
            delta{l} = zeros(N_neurons{l},1);
        end
        for test = 1:size_batch
            pattern = randi(length(training_set_out)); %test; %
            numbers(test) = pattern;
            V{1} = (training_set_in(pattern,:))';
            V_tot{1} = V_tot{1} + V{1}./size_batch;
            for l = 2:layers
                b_value = w{l}'*V{l-1} - theta{l};
                V{l} = g(b_value);
                V_tot{l} = V_tot{l} + V{l}./size_batch;
            end
            
            b_value = w{end}'*V{end-1} - theta{end};
            delta_temp{layers} = g_prim(b_value)*(training_set_out(pattern,:)-V{end})./size_batch;
            delta{layers} = delta{layers} + delta_temp{layers};
        end
        for l = layers:-1:3
            b_value = w{l-1}'*V_tot{l-2} - theta{l-1};
            delta_temp{l-1} = w{l}*delta{l}.*g_prim(b_value);
            delta{l-1} = delta{l-1} + delta_temp{l-1};
        end
        for l = 2:layers
            delta_w = mu*V_tot{l-1}*delta{l}';
            w{l} = w{l} + delta_w;
            
            delta_theta = -mu*delta{l};
            theta{l} = theta{l} + delta_theta;
        end
    end
    
    
    % VALIDATION
    for pattern = 1:size(validation_set_in,1)
        V{1} = (validation_set_in(pattern,:))';
        for l = 2:layers
            b_value = w{l}'*V{l-1} - theta{l};
            V{l} = g(b_value);
        end
        O(pattern,:) = V{end};
        if O(pattern,:) == 0
            O = 1;
        end
    end
    C = (1/(2*length(validation_set_out)))*sum(abs(sign(O)-validation_set_out));
    if C<C_best
        C_best = C
        
    end
    if C<0.12
        w_best = w;
        theta_best = theta;
        break
    end
    if sum(ismember(C_vector, 0.1516)) == 7
        break
    end
end

find(abs(sign(O)-validation_set_out)==2)



