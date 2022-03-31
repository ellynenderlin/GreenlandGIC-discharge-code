function error = testmodelerror(modelfitobj, testdata)
% Evaluate a model's error (residual) using cross-validation testing data
% INPUT:    modelfitobj = the fit object from MATLAB fit function
%           testdata = matrix with x-values in column 1 and y-values in column 2
% OUTPUTS:  error = average residual
% Syntax: error = testmodelerror(modelfitobj, traindata)
y=testdata(:,2); % grav the actual y-values
ymod=feval(modelfitobj, testdata(:,1)); % evaluate model at the training data x-values
error=mean(abs(y-ymod)); % calculate RMSE
