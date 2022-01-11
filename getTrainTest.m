function [trainset, testset] = getTrainTest(data,pTrain)
% partition dataset in training and test sets
% INPUT: data = [ns,nd] dataset with ns samples, nd dimensions
%        pTrain = fraction of data for training 0-1
% OUTPUTS: trainset = data to use to train the model
%          testset = data for testing model
% Syntax: [trainset, testset] = getTrainTest(data,pTrain)
[ns,nd]=size(data); % ns=number of samples, nd= number of dimensons
nsTrain=round(pTrain*ns); % number of samples to train
Ix=1:ns; % vector od indicies to observations
IxTrain=randsample(Ix,nsTrain); % randomly grab rows in data
Ix3=ones(size(Ix)); % vector of ones
Ix3(IxTrain)=0; % set rows from training data to 0
IxTest=find(Ix3); % indecies to the test observations 
trainset=data(IxTrain,:); % data for training
testset=data(IxTest,:); % data for testing