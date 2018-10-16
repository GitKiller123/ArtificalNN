clear all; clc;

training_set = load('training_set.csv');
training_set_in = training_set(:,1:2);
training_set_out = training_set(:,3);

validation_set = load('validation_set.csv');
validation_set_in = validation_set(:,1:2);
validation_set_out = validation_set(:,3);

T = 10^5;
eta = 0.02;
size_batch = 200;
M1 = 3; %Amount of neurons in 1st hidden layer
M2 = 2; %Amount of neurons in 2nd hidden layer

g = @(b)(tanh(b));  %Activation function
g_prim = @(b)(sech(b)^2); %Derivative of the activation function
b = @(w,V,theta,i,l)(w{l}(:,i)'*V{l-1}' - theta{l}(i));
delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l-1)));

N_neurons = {2, M1, M2, 1}; %Amount of Neurons in each layer

layers = length(N_neurons); %Amount of layers
% for t = 1:100
for l = 2:length(N_neurons)
    w{l} = normrnd(0,1,N_neurons{l-1},N_neurons{l}); %Creating starting weights for each layer
    theta{l} = 0; %Creating starting thresholds for each layer
end
%
H_best = 10000;
C_best = 10000;
% for pattern = 1:size(validation_set_in,1)
%     V{1} = validation_set_in(pattern,:);
%     for l = 2:layers
%         for i = 1:N_neurons{l}
%             V{l}(i) = g(b(w,V,theta,i,l));
%         end
%     end
%     O(pattern,:) = V{end}(:);
% end
% C = (1/(2*length(validation_set_out)))*sum(abs(sign(O)-validation_set_out))
% end
for trial = 1:100 %size(training_set_in,1)
    % TRAINING
%     patterns = randi(size(training_set,1),80,1);
%     for p_counter = 1:length(patterns)
%         pattern = patterns(p_counter);
        pattern_num = randi(size(training_set,1)-size_batch);
        for pattern = pattern_num:pattern_num+size_batch
            p_counter = pattern;
        V{1} = training_set_in(pattern,:);
        for l = 2:layers
            for i = 1:N_neurons{l}
                V{l}(i) = g(b(w,V,theta,i,l));
            end
        end

        for i = 1:N_neurons{end}
            delta{layers}(i) = g_prim(b(w,V,theta,i,l))*(training_set_out(pattern,i)-V{end}(i));
        end
        for l = layers:-1:3
            for j = 1:N_neurons{l-1}
                for i = 1:N_neurons{l}
                    delta_temp(i) = delta_func(w,V,theta,delta,i,j,l);
                    %                 delta{l-1}(j) = delta{l}(i)*w{l}(j)*g_prim(b(w,V,theta,j,l));
% delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l)));
% b = @(w,V,theta,i,l)(0.5*(w{l}(:,i)'*V{l-1}' - theta{l}(i)));

                end
                delta{l-1}(j) = sum(delta_temp);
                delta_temp = [];
            end
        end
        for l = 2:layers-1
            delta_w = eta*delta{l}'*V{l-1};
%             delta_w = eta*delta{l}'*V{l};
            w{l} = w{l} + delta_w';
            
            delta_theta = -eta*delta{l}';
            theta{l} = theta{l} + delta_theta;
        end
        O_training(p_counter) = V{end};
        curr_training(p_counter) = training_set_out(pattern);
    end
    H = 0.5*(curr_training-O_training)*(curr_training-O_training)'/size_batch;

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
end

find(abs(sign(O)-validation_set_out)==2)


