% k-nearest neighbours
% 
% Pavel Trutman
% INRIA, 2016

deepagpaths;

%val_file='../../data/paris/features_F_nsynthrep.mat';
%tr_file='../../data/paris/features_F_nsynthrep_sample_1K.mat';

% prepare data
clear Xtr Ytr Ntr Xval Yval Nval;
data = matfile(val_file);
Xtr = data.tr_coefs;
Ytr = data.tr_f;
Ntr = data.tr_norm;
data_sample = matfile(tr_file);
Xval = data_sample.val_coefs;
Yval = data_sample.val_f;
Nval = data_sample.val_norm;

% properties
k = 1;
noise_std = 1/1000;
noise=true;

% data normalize
X_mean = mean(Xtr, 2);
X_std = std(Xtr, [], 2);
std_nz = find(X_std ~= 0);
Xtr = Xtr(std_nz, :);
Xval = Xval(std_nz, :);
X_mean = X_mean(std_nz);
X_std = X_std(std_nz);
Xtr = (Xtr - repmat(X_mean, 1, size(Xtr, 2)))./repmat(X_std, 1, size(Xtr, 2));
Xval = (Xval - repmat(X_mean, 1, size(Xval, 2)))./repmat(X_std, 1, size(Xval, 2));

%{
% uncomment this to remove linear dependent elements, it helps when the
% covariance matrix is badly conditioned

C = cov(Xtr');
keep = [];
removed = [];
while size(keep, 1) + size(removed, 1) ~= size(C, 1)
  [~, idxs] = sort(sum(abs(C), 2), 'descend');
  idxs = setdiff(setdiff(idxs, removed, 'stable'), keep, 'stable');
  fprintf([num2str(rcond(C([idxs; keep], [idxs; keep]))), '\n']);
  idx = idxs(1);
  if any(abs(C(idx, setdiff(idxs, idx))) > 0.999)
    removed = [removed; idx];
  else
    keep = [keep; idx];
  end
end
fprintf([num2str(rcond(C(keep, keep))), '\n']);
Xtr = Xtr(keep, :);
Xval = Xval(keep, :);
C = cov(Xtr');
Cinv = inv(C);
%}

% noise
if noise ~= 0
  noise = normrnd(0, noise_std, size(Xval));
  Xval = Xval + noise;
end

val_error = zeros(1, size(Xval, 2));
val_error_k = zeros(k, size(Xval, 2));
for i = 1:size(Xval, 2)
  repS = adprintf({}, [num2str(i), '/', num2str(size(Xval, 2))]);
  YvalDN = Yval(:, i).*Nval(:, i);
  dist = Xtr - repmat(Xval(:, i), 1, size(Xtr, 2));
  dist = sum(dist.^2, 1); % Euclidean distance
  %dist = sum(dist.*(Cinv*dist), 1); % Mahalanobis distance  %#ok<MINV>
  Yoval = zeros(2, k);
  for kk = 1:k
    [~, idx] = min(dist);
    dist(idx) = Inf;
    Yoval(:, kk) = Ytr(:, idx).*Ntr(:, idx);
    val_error_k(kk, i) = norm(YvalDN - Yoval(:, kk));
  end
  Yoval = mean(Yoval, 2);
  val_error(i) = norm(YvalDN - Yoval);
  rmprintf(repS);
end

% save result
[pathstr, name, ~] = fileparts(mfilename('fullpath'));
save([pathstr, '/results/', name, '-' datestr(now, 'yymmdd-HHMMSS') '.mat'], 'k', 'val_error', 'val_error_k');
clear X_mean X_std Xtr Xval Ytr Yval YvalDN Ntr Nval val_error val_error_k repS dist k kk i idx Yoval name pathstr;