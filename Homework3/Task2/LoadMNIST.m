function [xTrain, tTrain, xValid, tValid, xTest, tTest] = LoadMNIST(exerciseNumber)

% This MATLAB function is an adaptation of the function prepareData.m,
% copyrighted in 2018 by The MathWorks, Inc. It's intended use is for the
% course Artificial Neural Networks (FFR135/FIM720) at Chalmers University
% of Technology/University of Gothenburg, Gothenburg, Sweden, 2018.
%
%%%%%%%%%
% Input %
%%%%%%%%%
%
% exerciseNumber: A number that determines the format of the returned data.
% This number may be 1, 2, 3 or 4. See separate course material for a
% documentation.
%
%%%%%%%%%%
% Output %
%%%%%%%%%%
% 
% xTrain: The input patterns of the training set
% tTrain: The output patterns of the training set
% xValid: The input patterns of the validation set
% tValid: The output patterns of the validation set
% xTest: The input patterns of the test set
% tTest: The output patterns of the test set
% 
% The format of the output data depends on the in-argument provided as
% "exerciseNumber". See separate course material for a documentation.

if ismember(exerciseNumber,[1 2])
    array2DFormat = true;
elseif ismember(exerciseNumber,[3 4])
    array2DFormat = false;
else
    assert(false,'"Exercise number" must equal 1, 2, 3 or 4.');
end

% Check for the existence of the MNIST files and download them if necessary

mnistDir = 'mnist_data';
files = {   "train-images-idx3-ubyte",...
            "train-labels-idx1-ubyte",...
            "t10k-images-idx3-ubyte",...
            "t10k-labels-idx1-ubyte"  };

% boolean for testing if the files exist
% basically, check for existence of "data" directory
download = exist(fullfile(pwd, mnistDir), 'dir') ~= 7;

if download
    disp('Downloading files...')
    mkdir(mnistDir)
    webPrefix = "http://yann.lecun.com/exdb/mnist/";
    webSuffix = ".gz";

    filenames = files + webSuffix;
    for ii = 1:numel(files)
        websave(fullfile(mnistDir, filenames{ii}),...
            char(webPrefix + filenames(ii)));
    end
    disp('Download complete.')
    
    % unzip the files
    cd(mnistDir)
    gunzip *.gz
    
    % return to main directory
    cd ..
end

% Extract the MNIST images into arrays

disp('Preparing MNIST data...');

% Read headers for training set image file
fid = fopen(fullfile(mnistDir, char(files{1})), 'r', 'b');
magicNum = fread(fid, 1, 'uint32');
numImgs  = fread(fid, 1, 'uint32');
numRows  = fread(fid, 1, 'uint32');
numCols  = fread(fid, 1, 'uint32');

% Read the data part 
rawImgDataTrain = uint8(fread(fid, numImgs * numRows * numCols, 'uint8'));
fclose(fid);

% Reshape the data part into a 4D array
rawImgDataTrain = reshape(rawImgDataTrain, [numRows, numCols, numImgs]);
rawImgDataTrain = permute(rawImgDataTrain, [2,1,3]);
imgDataTrain(:,:,1,:) = uint8(rawImgDataTrain(:,:,:));

% Read headers for training set label file
fid = fopen(fullfile(mnistDir, char(files{2})), 'r', 'b');
magicNum  = fread(fid, 1, 'uint32');
numLabels = fread(fid, 1, 'uint32');

% Read the data for the labels
labelsTrain = fread(fid, numLabels, 'uint8');
fclose(fid);

% Process the labels
labelsTrain = categorical(labelsTrain);

% Read headers for test set image file
fid = fopen(fullfile(mnistDir, char(files{3})), 'r', 'b');
magicNum = fread(fid, 1, 'uint32');
numImgs  = fread(fid, 1, 'uint32');
numRows  = fread(fid, 1, 'uint32');
numCols  = fread(fid, 1, 'uint32');

% Read the data part 
rawImgDataTest = uint8(fread(fid, numImgs * numRows * numCols, 'uint8'));
fclose(fid);

% Reprocess the data part into a 4D array
rawImgDataTest = reshape(rawImgDataTest, [numRows, numCols, numImgs]);
rawImgDataTest = permute(rawImgDataTest, [2,1,3]);
imgDataTest = uint8(zeros(numRows, numCols, 1, numImgs));
imgDataTest(:,:,1,:) = uint8(rawImgDataTest(:,:,:));

% Read headers for test set label file
fid = fopen(fullfile(mnistDir, char(files{4})), 'r', 'b');
magicNum  = fread(fid, 1, 'uint32');
numLabels = fread(fid, 1, 'uint32');

% Read the data for the labels
labelsTest = fread(fid, numLabels, 'uint8');
fclose(fid);

% Process the labels
labelsTest = categorical(labelsTest);

disp('MNIST data preparation complete.');

% Split the MNIST training set into a validation set and a proper training set

rngState = rng; % store this state of the random number generator
rng(123); % set the random number generator
idxValidationSet = randperm(60000,10000);
isValidationSet = ismember(1:60000,idxValidationSet);
rng(rngState); % reset the random number generator to the original state

% Up to now, we have loaded precisely the data as the original
% prepareData.m function. The variables are:
% imgDataTrain, labelsTrain, imgDataTest, labelsTest

% Create our output sets
xTrain = imgDataTrain(:,:,:,not(isValidationSet));
tTrain = labelsTrain(not(isValidationSet));
xValid = imgDataTrain(:,:,:,isValidationSet);
tValid = labelsTrain(isValidationSet);
xTest = imgDataTest;
tTest = labelsTest;

if ~array2DFormat
    return;
end

% Change format of output sets

xTrain = double(xTrain);
tTrain = double(tTrain);
tmp = zeros(10,50000);
for i = 1:50000
    tmp(tTrain(i),i) = 1;
end
tTrain = tmp;
xTrain = reshape(xTrain,[784 50000]);
xTrain = xTrain/255;

xValid = double(xValid);
tValid = double(tValid);
tmp = zeros(10,10000);
for i = 1:10000
    tmp(tValid(i),i) = 1;
end
tValid = tmp;
xValid = reshape(xValid,[784 10000]);
xValid = xValid/255;

xTest = double(xTest);
tTest = double(tTest);
tmp = zeros(10,10000);
for i = 1:10000
    tmp(tTest(i),i) = 1;
end
tTest = tmp;
xTest = reshape(xTest,[784 10000]);
xTest = xTest/255;


end % end of function
