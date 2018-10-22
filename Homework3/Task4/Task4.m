clear all; clc;
tic

[xTrain, tTrain, xValid, tValid, xTest, tTest] = LoadMNIST(4);
% imshow(mat2gray(xTrain(:,:,1,25081)))
curr_network = 2;


switch curr_network
    case 1
        layers = [imageInputLayer([28 28 1])
            convolution2dLayer(5,20,'Padding',1,'Stride',1)
            reluLayer
            maxPooling2dLayer(2,'Stride',2,'Padding',0)
            fullyConnectedLayer(100)
            reluLayer
            fullyConnectedLayer(10)
            softmaxLayer
            classificationLayer];
        options = trainingOptions('sgdm',...
            'MiniBatchSize', 8192,...
            'ValidationData', {xValid, tValid},...
            'ValidationFrequency', 30,...
            'MaxEpochs',60,...
            'Plots', 'Training-Progress',...
            'L2Regularization', 0, ...
            'Momentum', 0.9, ...
            'ValidationPatience', 5, ...
            'Shuffle', 'every-epoch', ...
            'InitialLearnRate', 0.001);
        [net info] = trainNetwork(xTrain, tTrain, layers, options);
    case 2
        layers = [imageInputLayer([28 28 1])
            convolution2dLayer(3,20,'Padding',1,'Stride',1)
            batchNormalizationLayer
            reluLayer
            maxPooling2dLayer(2,'Stride',2,'Padding',0)
            convolution2dLayer(3,30,'Padding',1,'Stride',1)
            batchNormalizationLayer
            reluLayer
            maxPooling2dLayer(2,'Stride',2,'Padding',0)
            convolution2dLayer(3,50,'Padding',1,'Stride',1)
            batchNormalizationLayer
            reluLayer
            fullyConnectedLayer(10)
            softmaxLayer
            classificationLayer];
        options = trainingOptions('sgdm',...
            'MiniBatchSize', 8192,...
            'ValidationData', {xValid, tValid},...
            'ValidationFrequency', 30,...
            'MaxEpochs',30,...
            'Plots', 'Training-Progress',...
            'L2Regularization', 0, ...
            'Momentum', 0.9, ...
            'ValidationPatience', 5, ...
            'Shuffle', 'every-epoch', ...
            'InitialLearnRate', 0.01);
        [net info] = trainNetwork(xTrain, tTrain, layers, options);
end
%% Classification error

O_train = net.classify(xTrain);
error_count = 0;
for i = 1:size(xTrain,4)
    if tTrain(i) ~= O_train(i)
        error_count = error_count + 1;
    end
end
C_train = (1/size(xTrain,4))*error_count;
O_valid = net.classify(xValid);
error_count = 0;
for i = 1:size(xValid,4)
    if tValid(i) ~= O_valid(i)
        error_count = error_count + 1;
    end
end
C_valid = (1/size(xValid,4))*error_count;
O_test = net.classify(xTest);
error_count = 0;
for i = 1:size(xTest,4)
    if tTest(i) ~= O_test(i)
        error_count = error_count + 1;
    end
end
C_test = (1/size(xValid,4))*error_count;

