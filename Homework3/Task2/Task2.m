clear all; clc;
tic

[xTrain, tTrain, xValid, tValid, xTest, tTest] = LoadMNIST(2);
% imshow(reshape(xTrain(:,25081),[28 28]));

training_set_in = xTrain';
training_set_out = tTrain';

validation_set_in = xValid';
validation_set_out = tValid';

test_set_in = xTest';
test_set_out = tTest';

total_epochs = 50;
mu = 3e-3;
epoch_size = size(training_set_in,1);
size_batch = 10;

g = @(b)((1+exp(-b)).^(-1));  %Activation function
g_prim = @(b)(g(b)-g(b).^2); %Derivative of the activation function
x_avg = mean(training_set_in,1)';

M = 30; %Number of neurons
w = {};
theta = {};
N_neurons = {784, M, M, M, M, 10}; %Amount of Neurons in each layer
layers = length(N_neurons); %Amount of layers
for l = 2:length(N_neurons)
    for i = 1:N_neurons{l}
        w{l}(:,i) = normrnd(0, 1/(sqrt(N_neurons{l-1})),[N_neurons{l-1} 1]); %Creating starting weights for each layer
    end
    theta{l} = zeros(N_neurons{l},1); %Creating starting thresholds for each layer
end
counter2 = 1;
for epoch = 1:total_epochs
    if mod(epoch,5) == 0
        epoch
    end
    sequence = randperm(epoch_size);
    counter = 1;
    
    for trial = 1:epoch_size/size_batch
        for l = 2:layers
            delta_w{l} = zeros(N_neurons{l-1},N_neurons{l});
            delta_theta{l} = zeros(N_neurons{l},1);
            delta_tot{l} = zeros(N_neurons{l},1);
        end
        
        dH = 0;
        for test = 1:size_batch
            pattern = sequence(counter);
            numbers(test) = pattern;
            V{1} = (training_set_in(pattern,:))'-x_avg;
            for l = 2:layers
                b_value = w{l}'*V{l-1} - theta{l};
                V{l} = g(b_value);
            end
            b_value = w{end}'*V{end-1} - theta{end};
            delta{layers} = g_prim(b_value).*(training_set_out(pattern,:)'-V{end});
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
        for l = 2:layers
        u_l_tot{l} = zeros(N_neurons{l},1);
    end
    for pattern = 1:size(training_set_in,1)
        V{1} = (training_set_in(pattern,:))'-x_avg;
        for l = 2:layers
            b_value = w{l}'*V{l-1} - theta{l};
            V{l} = g(b_value);
        end
        b_value = w{end}'*V{end-1} - theta{end};
        u_l{layers} = g_prim(b_value).*(training_set_out(pattern,:)'-V{end});
        u_l_tot{layers}(pattern) = norm(u_l{layers});
        for l = layers:-1:3
            b_value = w{l-1}'*V{l-2} - theta{l-1};
            u_l{l-1} = w{l}*u_l{l}.*g_prim(b_value);
            u_l_tot{l-1}(pattern) = norm(u_l{l-1});
        end
        pos = find(V{end}==max(V{end}));
        V{end}(pos) = 1;
        V{end}([1:pos-1 pos+1:end]) = 0;
        H_temp(pattern) = (tTrain(:,pattern)-V{end})'*(tTrain(:,pattern)-V{end});
    end
    for l = 2:layers
        U{l}(epoch) = sum(u_l_tot{l});
    end
    H(epoch) = 0.5*sum(H_temp);
end

toc
%% plot
figure(1)
for l = 2:layers
    semilogy(U{l})
    hold on
    xlabel('Epochs')
    ylabel('Learning speed')
end
legend('layer 2', 'layer 3', 'layer 4', 'layer 5', 'layer 6')
figure(2)
plot(H)
