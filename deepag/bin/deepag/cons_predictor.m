% Constant predictor
% 
% Pavel Trutman
% INRIA, 2016

% prepare data
data = matfile('../../data/paris/features.mat');
Ytr = data.tr_f;
Ntr = data.tr_norm;
data = matfile('../../data/paris/features_sample_100k.mat');
Yval = data.val_f;
Nval = data.val_norm;
Ytst = data.tst_f;
Ntst = data.tst_norm;
clear data;

% compute prediction
Ypred = mean(Ytr.*Ntr, 2);

% compute error
tr_error  = mean(sqrt(sum((Ytr.*Ntr   - repmat(Ypred, 1, size(Ytr,  2))).^2)));
val_error = mean(sqrt(sum((Yval.*Nval - repmat(Ypred, 1, size(Yval, 2))).^2)));
tst_error = mean(sqrt(sum((Ytst.*Ntst - repmat(Ypred, 1, size(Ytst, 2))).^2)));