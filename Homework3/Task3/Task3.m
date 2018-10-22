clear all; clc;
tic

[xTrain, tTrain, xValid, tValid, xTest, tTest] = LoadMNIST(3);
% imshow(mat2gray(xTrain(:,:,1,25081)))
total_networks = 3;

for curr_network = 1:total_networks
    if curr_network == 1
        layers = [imageInputLayer([28 28 1])
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(10)
            softmaxLayer
            classificationLayer];
        options = trainingOptions('sgdm',...
            'MiniBatchSize', 8192,...
            'ValidationData', {xValid, tValid},...
            'ValidationFrequency', 30,...
            'MaxEpochs',200,...
            'Plots', 'Training-Progress',...
            'L2Regularization', 0, ...
            'Momentum', 0.9, ...
            'ValidationPatience', 5, ...
            'Shuffle', 'every-epoch', ...
            'InitialLearnRate', 0.01);
        [net{curr_network} info{curr_network}] = trainNetwork(xTrain, tTrain, layers, options);
    elseif curr_network == 2
        layers = [imageInputLayer([28 28 1])
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(10)
            softmaxLayer
            classificationLayer];
        options = trainingOptions('sgdm',...
            'MiniBatchSize', 8192,...
            'ValidationData', {xValid, tValid},...
            'ValidationFrequency', 30,...
            'MaxEpochs',200,...
            'Plots', 'Training-Progress',...
            'L2Regularization', 0, ...
            'Momentum', 0.9, ...
            'ValidationPatience', 5, ...
            'Shuffle', 'every-epoch', ...
            'InitialLearnRate', 0.01);
        [net{curr_network} info{curr_network}] = trainNetwork(xTrain, tTrain, layers, options);
    elseif curr_network == 3
        layers = [imageInputLayer([28 28 1])
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(10)
            softmaxLayer
            classificationLayer];
        options = trainingOptions('sgdm',...
            'MiniBatchSize', 8192,...
            'ValidationData', {xValid, tValid},...
            'ValidationFrequency', 30,...
            'MaxEpochs',200,...
            'Plots', 'Training-Progress',...
            'L2Regularization', 0.03, ...
            'Momentum', 0.9, ...
            'ValidationPatience', 5, ...
            'Shuffle', 'every-epoch', ...
            'InitialLearnRate', 0.01);
        [net{curr_network} info{curr_network}] = trainNetwork(xTrain, tTrain, layers, options);
    end
    %% Classification error
    
    O_train = net{curr_network}.classify(xTrain);
    error_count = 0;
    for i = 1:size(xTrain,4)
        if tTrain(i) ~= O_train(i)
            error_count = error_count + 1;
        end
    end
    C_train(curr_network) = (1/size(xTrain,4))*error_count;
    O_valid = net{curr_network}.classify(xValid);
    error_count = 0;
    for i = 1:size(xValid,4)
        if tValid(i) ~= O_valid(i)
            error_count = error_count + 1;
        end
    end
    C_valid(curr_network) = (1/size(xValid,4))*error_count;
    O_test = net{curr_network}.classify(xTest);
    error_count = 0;
    for i = 1:size(xTest,4)
        if tTest(i) ~= O_test(i)
            error_count = error_count + 1;
        end
    end
    C_test(curr_network) = (1/size(xValid,4))*error_count;
end

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
%         sequence = randperm(epoch_size);
%         counter = 1;
%         for trial = 1:epoch_size/size_batch
%             for l = 2:layers
%                 delta_w{l} = zeros(N_neurons{l-1},N_neurons{l});
%                 delta_theta{l} = zeros(N_neurons{l},1);
%             end
%             dH = 0;
%             for test = 1:size_batch
%                 pattern = sequence(counter);
%                 numbers(test) = pattern;
%                 V{1} = (training_set_in(pattern,:))'-x_avg;
%                 for l = 2:layers
%                     b_value = w{l}'*V{l-1} - theta{l};
%                     V{l} = g(b_value);
%                 end
%                 O_H(test,:) = V{end};
%                 b_value = w{end}'*V{end-1} - theta{end};
%                 delta{layers} = g_prim(b_value).*(training_set_out(pattern,:)'-V{end});
%                 for l = layers:-1:3
%                     b_value = w{l-1}'*V{l-2} - theta{l-1};
%                     delta{l-1} = w{l}*delta{l}.*g_prim(b_value);
%                 end
%                 for l = 2:layers
%                     delta_w{l} = delta_w{l} + mu*V{l-1}*delta{l}';
%                     delta_theta{l} = delta_theta{l} - mu*delta{l};
%                 end
%                 counter = counter + 1;
%             end
%             for l = 2:layers
%                 w{l} = w{l} + delta_w{l};
%                 theta{l} = theta{l} + delta_theta{l};
%             end
%         end
%     end
%     C_train_tot{curr_case} = C_train;
%     C_val_tot{curr_case} = C_val;
% end
% %% Test set
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
% toc
% %% plot
% for i = 1:total_cases
%     subplot(2,2,i)
%     semilogy(C_train_tot{i},'b')
%     hold on
%     semilogy(C_val_tot{i},'-r')
%     axis([0 30 0.001 1])
%     xlabel('Epochs')
%     ylabel('Classification Error')
%     legend('Training set', 'Validation set')
%     title('Network ' + string(i))
% end
