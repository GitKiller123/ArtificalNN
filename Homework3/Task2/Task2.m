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
% b = @(w,V,theta,i,l)(w{l}(:,i)'*V{l-1} - theta{l}(i));
% delta_func = @(w,V,theta,delta,i,j,l)(delta{l}(i)*w{l}(j,i)*g_prim(b(w,V,theta,j,l-1)));

total_cases = 4;
x_avg = mean(training_set_in,1)';

% for curr_case = 1:total_cases
%     curr_case
M = 30; %Number of neurons
    w = {};
    theta = {};
    N_neurons = {784, M, M, M, M, 10}; %Amount of Neurons in each layer
    layers = length(N_neurons); %Amount of layers
%     C_best.val(curr_case) = 1;
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
%         for pattern = 1:size(training_set_in,1)
%             V{1} = (training_set_in(pattern,:))'-x_avg;
%             for l = 2:layers
%                 b_value = w{l}'*V{l-1} - theta{l};
%                 V{l} = g(b_value);
%             end
%             pos = find(V{end}==max(V{end}));
%             V{end}(pos) = 1;
%             V{end}([1:pos-1 pos+1:end]) = 0;
%             O_train(pattern,:) = V{end};
%         end
%         C_train(epoch) = (1/(2*size(training_set_out,1)))*sum(sum(abs(training_set_out-O_train)));
%         % VALIDATION
%         for pattern = 1:size(validation_set_in,1)
%             V{1} = (validation_set_in(pattern,:))'-x_avg;
%             for l = 2:layers
%                 b_value = w{l}'*V{l-1} - theta{l};
%                 V{l} = g(b_value);
%             end
%             pos = find(V{end}==max(V{end}));
%             V{end}(pos) = 1;
%             V{end}([1:pos-1 pos+1:end]) = 0;
%             O_val(pattern,:) = V{end};
%         end
%         C_val(epoch) = (1/(2*size(validation_set_out,1)))*sum(sum(abs(validation_set_out-O_val)));
%         if C_val(epoch)<C_best.val(curr_case)
%             C_best.val(curr_case) = C_val(epoch);
%             C_best.train(curr_case) = C_train(epoch);
%             C_best.epoch(curr_case) = epoch;
%             C_best.(['weights_' num2str(curr_case)]) = w;
%             C_best.(['theta_' num2str(curr_case)]) = theta;
%         end
        sequence = randperm(epoch_size);
        counter = 1;
        for trial = 1:epoch_size/size_batch
            for l = 2:layers
                delta_w{l} = zeros(N_neurons{l-1},N_neurons{l});
                delta_theta{l} = zeros(N_neurons{l},1);
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
                O_H(test,:) = V{end};
                b_value = w{end}'*V{end-1} - theta{end};
                delta{layers} = g_prim(b_value).*(training_set_out(pattern,:)'-V{end});
                for l = layers:-1:3
                    b_value = w{l-1}'*V{l-2} - theta{l-1};
                    delta{l-1} = w{l}*delta{l}.*g_prim(b_value);
                end
                for l = 2:layers
                    delta_w{l} = delta_w{l} + mu*V{l-1}*delta{l}';
                    delta_theta{l} = delta_theta{l} - mu*delta{l};
                    U{l}(counter2) = norm(delta{l});
                end
                
                counter = counter + 1;
            end
            for l = 2:layers
                w{l} = w{l} + delta_w{l};
                theta{l} = theta{l} + delta_theta{l};
            end
        end
        for l = 2:layers
            U{l}(counter2) = norm(delta{l});
        end
        counter2 = counter2 + 1;
    end
    C_train_tot{curr_case} = C_train;
    C_val_tot{curr_case} = C_val;
% end
%% Test set
% for curr_case = 1:total_cases
%     w = C_best.(['weights_' num2str(curr_case)]);
%     theta = C_best.(['theta_' num2str(curr_case)]);
%     layers = length(w);
%     for pattern = 1:size(test_set_in,1)
%         V{1} = (test_set_in(pattern,:))'-x_avg;
%         for l = 2:layers
%             b_value = w{l}'*V{l-1} - theta{l};
%             V{l} = g(b_value);
%         end
%         pos = find(V{end}==max(V{end}));
%         V{end}(pos) = 1;
%         V{end}([1:pos-1 pos+1:end]) = 0;
%         O_test(pattern,:) = V{end};
%     end
%     C_best.test(curr_case) = (1/(2*size(test_set_out,1)))*sum(sum(abs(test_set_out-O_test)));
% end
% save('Task1_1_variables','C_best','C_val','C_train')
toc
%% plot
for l = 2:layers
    plot(U{l})
    hold on
    xlabel('Epochs')
    ylabel('Learning speed')
end
legend('layer 2', 'layer 3', 'layer 4', 'layer 5', 'layer 6')
