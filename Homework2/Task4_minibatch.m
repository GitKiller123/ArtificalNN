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
M1 = 4; %Amount of neurons in 1st hidden layer
M2 = 9; %Amount of neurons in 2nd hidden layer

g = @(b)(tanh(b));  %Activation function
g_prim = @(b)(sech(b)^2); %Derivative of the activation function
b = @(w,V,theta,i,l)(0.5*(w{l}(:,i)'*V{l-1}' - theta{l}(i)));
delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l-1)));

N_neurons = {2, M1, M2, 1}; %Amount of Neurons in each layer

layers = length(N_neurons); %Amount of layers

H_best = 10000;
C_best = 0.17;
C_vector = [];

for attempt = 1:20
    C_best = 0.17;
    C_vector = [];
    if attempt > 5
        M1 = M1 + 1; %Amount of neurons in 1st hidden layer
        M2 = M2 + 1; %Amount of neurons in 2nd hidden layer
        N_neurons = {2, M1, M2, 1}; %Amount of Neurons in each layer
    end
    for l = 2:length(N_neurons)
        w{l} = normrnd(0,0.2,N_neurons{l-1},N_neurons{l}); %Creating starting weights for each layer
        theta{l} = zeros(N_neurons{l},1); %Creating starting thresholds for each layer
    end
    for trial = 1:size_batch:size(training_set_in,1)
        % TRAINING
        V{1} = 0;
        if trial + size_batch > size(training_set_in,1)
            break
        end
        for l = 2:layers
            V{l} = [];
        end
%         jump = randi(size(training_set_in,1)-size_batch);
        for pattern = trial:trial+size_batch
            V_temp{1} = training_set_in(pattern,:);
            V{1} = V{1} + V_temp{1}/size_batch;
            for l = 2:layers
                for i = 1:N_neurons{l}
                    V_temp{l}(i) = g(b(w,V_temp,theta,i,l));
                end
                V{l} = V_temp{l}./size_batch;
            end
        end
        
        for i = 1:N_neurons{end}
            delta{layers}(i) = g_prim(b(w,V,theta,i,l))*(training_set_out(pattern,i)-V{end}(i));
        end
        for l = layers:-1:3
            for j = 1:N_neurons{l-1}
                for i = 1:N_neurons{l}
                    delta_temp(i) = delta_func(w,V,theta,delta,i,j,l);                    
                end
                delta{l-1}(j) = sum(delta_temp);
                delta_temp = [];
            end
        end
        for l = 2:layers-1
            delta_w = eta*delta{l}'*V{l-1};
            w{l} = w{l} + delta_w';
            
            delta_theta = -eta*delta{l}';
            theta{l} = theta{l} + delta_theta;
        end
        
        % VALIDATION
        for pattern = 1:size(validation_set_in,1)
            V{1} = validation_set_in(pattern,:);
            for l = 2:layers
                for i = 1:N_neurons{l}
                    V{l}(i) = g(b(w,V,theta,i,l));
                end
            end
            O(pattern,:) = V{end}(:);
        end
        C = (1/(2*length(validation_set_out)))*sum(abs(sign(O)-validation_set_out))
        if C<C_best
            C_best = C
        end
        if C<0.12
            w_best = w;
            theta_best = theta;
            break
        end
        if C_vector<8
            C_vector = [C_vector C];
        else
            C_vector = C_vector(2:end);
            C_vector = [C_vector C];
        end
        if sum(ismember(C_vector, 0.1516)) == 7
            break
        end
    end
    if C<0.12
        w_best = w;
        theta_best = theta;
        break
    end
end

find(abs(sign(O)-validation_set_out)==2)


