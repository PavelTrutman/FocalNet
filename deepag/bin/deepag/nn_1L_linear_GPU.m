% 1-layer linear neural network
% 
% Pavel Trutman
% INRIA, 2016

% init matconvnet
matconvnet_setup_GPU;

% prepare data
clear Xtr Ytr Ntr Xval Yval Nval;
data = matfile('../../data/paris/features.mat');
Xtr(1, 1, :, :) = gpuArray(single(data.tr_coefs));
Ytr = gpuArray(data.tr_f);
Ntr = gpuArray(data.tr_norm);
data = matfile('../../data/paris/features_sample_100k.mat');
Xval(1, 1, :, :) = gpuArray(single(data.val_coefs));
Yval = gpuArray(data.val_f);
Nval = gpuArray(data.val_norm);
clear data;

% properties
di = size(Xtr, 3); % input dimension
do = size(Ytr, 1); % output dimension
lrate = 3e-8; % learning rate
mrate = 0.1; % momentum rate
nsamples = size(Xtr, 4);
nepoch = 30;
plotPeriod = nsamples/2;
batchSize = 100;

% data normalize
X_mean = mean(Xtr, 4);
X_std = std(Xtr, [], 4);
X_std(X_std == 0) = 1;
Xtr = (Xtr - repmat(X_mean, 1, 1, 1, size(Xtr, 4)))./repmat(X_std, 1, 1, 1, size(Xtr, 4));
Xval = (Xval - repmat(X_mean, 1, 1, 1, size(Xval, 4)))./repmat(X_std, 1, 1, 1, size(Xval, 4));
Y_mean = mean(Ytr, 2);
Y_std = std(Ytr, [], 2);
Y_std(Y_std == 0) = 1;

% randomly initialize parameters of the model
%w_inv = squeeze(Xtr)'\squeeze(Ytr)';
%w(1, 1, :, :) = single(w_inv);
w = single(normrnd(0, 1/sqrt(di), [1, 1, di, do]));
%w = single(zeros([1, 1, di, do]));
w_m = zeros(size(w));

clear net;
net.layers{1} = struct(...
    'name', 'linear1', ...
    'type', 'conv', ...
    'weights', {{w, []}}, ...
    'stride', 1, ...
    'pad', 0, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{2} = struct(...
    'name', 'denorm', ...
    'type', 'denorm', ...
    'norm', {{Y_mean, Y_std}}, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{3} = struct(...
    'name', 'loss', ...
    'type', 'lossL2', ...
    'norm', [], ...
    'target', [], ...
    'opts', {{}}, ...
    'precious', false);

tr_error = zeros(1, idivide(int32(nepoch*nsamples), plotPeriod, 'ceil') + 1);
val_error = zeros(1, idivide(int32(nepoch*nsamples), plotPeriod, 'ceil') + 1);

net = vl_simplenn_move(net, 'gpu');

j = 1;
tic;
for epoch = 1:nepoch
    
  % go randomly through samples
  r = randperm(nsamples);
  for batch = 1:(nsamples/batchSize);
      
    if mod(batch*batchSize, plotPeriod) == 0 || mod(batch, nsamples/batchSize) == 0
      
      % training error
      net.layers{end}.target = Ytr;
      net.layers{end}.norm = Ntr;
      res = vl_simplenn(net, Xtr, [], [], 'ConserveMemory', true, 'CuDNN', true);
      tr_error(j) = gather(res(end).x);
      clear res;
      
      % validation error
      net.layers{end}.target = Yval;
      net.layers{end}.norm = Nval;
      res = vl_simplenn(net, Xval, [], [], 'ConserveMemory', true, 'CuDNN', true);
      val_error(j) = gather(res(end).x);
      clear res;
      
      fprintf('%3i/%3i: log(tr): %7.3s, log(val): %7.3s, tr: %7.3s, val: %7.3s, sec: %5.2fs\n', epoch, nepoch, log10(tr_error(j)), log10(val_error(j)), tr_error(j), val_error(j), toc);
      tic;

      j = j + 1;
    end    
      
    n = r(((batch - 1)*batchSize + 1):(batch*batchSize));
    X = Xtr(:, :, :, n);
    Y = Ytr(:, n);
    N = Ntr(:, n);
    
    % Forward-Backward pass
    net.layers{end}.target = Y;
    net.layers{end}.norm = N;
    res = vl_simplenn(net, X, 1, [], 'ConserveMemory', true, 'CuDNN', true);
    
    % Update weights
    w_m = mrate*w_m + lrate*res(1).dzdw{1};
    net.layers{1}.weights{1} = net.layers{1}.weights{1} - w_m;
  end
end

% training error
net.layers{end}.target = Ytr;
net.layers{end}.norm = Ntr;
res = vl_simplenn(net, Xtr, [], [], 'ConserveMemory', true, 'CuDNN', true);
tr_error(j) = gather(res(end).x); %#ok<NASGU>
clear res;
      
% validation error
net.layers{end}.target = Yval;
net.layers{end}.norm = Nval;
res = vl_simplenn(net, Xval, [], [], 'ConserveMemory', true, 'CuDNN', true);
val_error(j) = gather(res(end).x); %#ok<NASGU>
clear res;

clear i j epoch n dsdy dsdw di do lrate mrate nepoch X Y Yhat Yotr Yoval;

plotPerEpoch = plotPeriod/nsamples; %#ok<NASGU>
[pathstr, name, ~] = fileparts(mfilename('fullpath'));
save([pathstr, '/results/', name, '-' datestr(now, 'yymmdd-HHMMSS') '.mat'], 'net', 'tr_error', 'val_error', 'plotPerEpoch');
clear tr_error val_error plotPerEpoch plotPeriod nsamples;
