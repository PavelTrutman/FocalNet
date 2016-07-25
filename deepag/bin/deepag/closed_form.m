% Pavel Trutman
% INRIA, 2016
% 
% Linear closed form solution

% prepare data
data = matfile('../../data/paris/features.mat');
Xtr = data.tr_coefs;
Ytr = data.tr_f;
Ntr = data.tr_norm;
data = matfile('../../data/paris/features_sample_100k.mat');
Xval = data.val_coefs;
Yval = data.val_f;
Nval = data.val_norm;
Xtst = data.tst_coefs;
Ytst = data.tst_f;
Ntst = data.tst_norm;
clear data;

W = Xtr'\Ytr';

YtrHat = transpose(Xtr'*W);
YvalHat = transpose(Xval'*W);
YtstHat = transpose(Xtst'*W);

tr_error  = mean(sqrt(sum((YtrHat.*Ntr   - Ytr.*Ntr  ).^2)));
val_error = mean(sqrt(sum((YvalHat.*Nval - Yval.*Nval).^2)));
tst_error = mean(sqrt(sum((YtstHat.*Ntst - Ytst.*Ntst).^2)));