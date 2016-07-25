% 2-layer neural network
% 
% Pavel Trutman
% INRIA, 2016

% init matconvnet
matconvnet_setup;

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
%lrate = 3e-8; % learning rate
lrate0 = 3e-8; % learning rate
mrate = 0.1; % momentum rate
nsamples = size(Xtr, 4);
nepoch = 20;
plotPeriod = nsamples/2;
batchSize = 100;
h1 = 500;
h2 = 500;

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
w1 = single(normrnd(0, 1/sqrt(di), [1, 1, di, h1]));
%w1 = single(zeros([1, 1, di, h1]));
w1_m = zeros(size(w1), 'gpuArray');
b1 = single(normrnd(0, 1/sqrt(di), [1, 1, h1]));
%b1 = single(zeros([1, 1, h1]));
b1_m = zeros(size(b1), 'gpuArray');
w2 = single(normrnd(0, 1/sqrt(h1), [1, 1, h1, h2]));
%w2 = single(zeros([1, 1, h1, h2]));
w2_m = zeros(size(w2), 'gpuArray');
b2 = single(normrnd(0, 1/sqrt(h1), [1, 1, h2]));
%b2 = single(zeros([1, 1, h2]));
b2_m = zeros(size(b2), 'gpuArray');
w3 = single(normrnd(0, 1/sqrt(h2), [1, 1, h2, do]));
%w3 = single(zeros([1, 1, h2, do]));
w3_m = zeros(size(w3), 'gpuArray');
b3 = single(normrnd(0, 1/sqrt(h2), [1, 1, do]));
%b3 = single(zeros([1, 1, do]));
b3_m = zeros(size(b3), 'gpuArray');

clear net;
net.layers{1} = struct(...
    'name', 'linear1', ...
    'type', 'conv', ...
    'weights', {{w1, b1}}, ...
    'stride', 1, ...
    'pad', 0, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{2} = struct(...
    'name', 'relu1', ...
    'type', 'relu', ...
    'leak', 0, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{3} = struct(...
    'name', 'linear2', ...
    'type', 'conv', ...
    'weights', {{w2, b2}}, ...
    'stride', 1, ...
    'pad', 0, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{4} = struct(...
    'name', 'relu2', ...
    'type', 'relu', ...
    'leak', 0, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{5} = struct(...
    'name', 'linear3', ...
    'type', 'conv', ...
    'weights', {{w3, b3}}, ...
    'stride', 1, ...
    'pad', 0, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{6} = struct(...
    'name', 'denorm', ...
    'type', 'denorm', ...
    'norm', {{Y_mean, Y_std}}, ...
    'opts', {{}}, ...
    'precious', false);
net.layers{7} = struct(...
    'name', 'loss', ...
    'type', 'lossL2', ...
    'target', [], ...
    'norm', [], ...
    'opts', {{}}, ...
    'precious', false);

tr_error = zeros(1, idivide(int32(nepoch*nsamples), plotPeriod, 'ceil'));
val_error = zeros(1, idivide(int32(nepoch*nsamples), plotPeriod, 'ceil'));

net = vl_simplenn_move(net, 'gpu');

fprintf('lrate0: %7.1s, mrate: %7.1s, batchSize: %3d\n', lrate0, mrate, batchSize);

j = 1;
tic;
for epoch = 1:nepoch
  
  lrate = lrate0/(1 + 4*epoch/nepoch);
    
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
      
      %{
      if j > 1          
        figure(1);
        semilogy((0:(j-1))*plotPeriod/nsamples, tr_error(1:j));
        hold on;
        semilogy((0:(j-1))*plotPeriod/nsamples, val_error(1:j));
        hold off;
        legend('Training dataset', 'Validating dataset', 'Location', 'southoutside');
        xlabel('Epoch');
        ylabel('Error [px]');
        xlim([0 (j-1)*plotPeriod/nsamples]);
        ylim([min([tr_error(1:j) val_error(1:j)]) max([tr_error(1:j) val_error(1:j)])]);
        grid on;
        drawnow;
      end
      %}
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
    % 1
    w1_m = mrate*w1_m + lrate*res(1).dzdw{1};
    net.layers{1}.weights{1} = net.layers{1}.weights{1} - w1_m;
    b1_m = mrate*b1_m + lrate*res(1).dzdw{2};
    net.layers{1}.weights{2} = net.layers{1}.weights{2} - b1_m;
    % 2
    w2_m = mrate*w2_m + lrate*res(3).dzdw{1};
    net.layers{3}.weights{1} = net.layers{3}.weights{1} - w2_m;
    b2_m = mrate*b2_m + lrate*res(3).dzdw{2};
    net.layers{3}.weights{2} = net.layers{3}.weights{2} - b2_m;
    % 3
    w3_m = mrate*w3_m + lrate*res(5).dzdw{1};
    net.layers{5}.weights{1} = net.layers{5}.weights{1} - w3_m;
    b3_m = mrate*b3_m + lrate*res(5).dzdw{2};
    net.layers{5}.weights{2} = net.layers{5}.weights{2} - b3_m;
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

clear i j r epoch batch n dsdy dsdw di do nepoch X Y Yhat Yotr Yoval;

plotPerEpoch = plotPeriod/nsamples; %#ok<NASGU>
[pathstr, name, ~] = fileparts(mfilename('fullpath'));
save([pathstr, '/results/', name, '-' datestr(now, 'yymmdd-HHMMSS') '.mat'], 'net', 'tr_error', 'val_error', 'plotPerEpoch', 'lrate', 'mrate', 'batchSize');
clear tr_error val_error  plotPerEpoch plotPeriod nsamples lrate mrate batchSize;