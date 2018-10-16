clear all; clc;

[xTrain, tTrain, xValid, tValid, xTest, tTest] = LoadMNIST(1);
% imshow(reshape(xTrain(:,25081),[28 28]));

training_set_in = xTrain';
training_set_out = tTrain';

validation_set_in = xValid';
validation_set_out = tValid';

x_avg = mean(xTrain,2);

T = 30;
mu = 0.03;
epoch = size(training_set_in,1);
size_batch = 10;
M1 = 784; %Amount of neurons in 1st hidden layer
M2 = 10; %Amount of neurons in 2nd hidden layer

g = @(b)((1+exp(b)).^(-1));  %Activation function
g_prim = @(b)(sech(b).^2); %Derivative of the activation function
% b = @(w,V,theta,i,l)(w{l}(:,i)'*V{l-1} - theta{l}(i));
delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l-1)));

N_neurons = {M1, M2}; %Amount of Neurons in each layer

layers = length(N_neurons); %Amount of layers

C_best = 1;

for l = 2:length(N_neurons)
    for i = 1:N_neurons{l}
        w{l}(:,i) = normrnd(0, 1/(sqrt(N_neurons{l-1})),[N_neurons{l-1} 1]); %Creating starting weights for each layer
    end
    theta{l} = zeros(N_neurons{l},1); %Creating starting thresholds for each layer
end
for epochs = 1:T
    sequence = randperm(epoch);
    counter = 1;
    for trial = 1:epoch/size_batch
        for l = 2:layers
            delta_w{l} = zeros(N_neurons{l-1},N_neurons{l});
            delta_theta{l} = zeros(N_neurons{l},1);
        end
        for test = 1:size_batch
            pattern = sequence(counter);
%             numbers(test) = pattern;
            V{1} = (training_set_in(pattern,:))';
            for l = 2:layers
                b_value = w{l}'*V{l-1} - theta{l};
                V{l} = g(b_value);
            end
            b_value = w{end}'*V{end-1} - theta{end};
            delta{layers} = g_prim(b_value).*(training_set_out(pattern,:)'-V{end})./size_batch;
            for l = layers:-1:3
                b_value = w{l-1}'*V{l-2} - theta{l-1};
                delta{l-1} = w{l}*delta{l}.*g_prim(b_value);
            end
            for l = 2:layers
                delta_w{l} = delta_w{l} + mu*V{l-1}*delta{l}';
                delta_theta{l} = delta_theta{l} - mu*delta{l};
            end
            counter = counter + 1;
        end
        for l = 2:layers
            w{l} = w{l} + delta_w{l};
            theta{l} = theta{l} + delta_theta{l};
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
    C = mean((1/(2*length(validation_set_out)))*sum(abs(round(O)-validation_set_out)))
    if C<C_best
        C_best = C
        
    end
    if C<0.12
        w_best = w;
        theta_best = theta;
        break
    end
%     if sum(ismember(C_vector, 0.1516)) == 7
%         break
%     end
end



