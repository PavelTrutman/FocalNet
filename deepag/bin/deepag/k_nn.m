% Pavel Trutman
% INRIA, 2016
% 
% k-nearest neighbours

deepagpaths;

% prepare data
clear Xtr Ytr Ntr Xval Yval Nval;
data = matfile('../../data/paris/features.mat');
Xtr = data.tr_coefs;
Ytr = data.tr_f;
Ntr = data.tr_norm;
data_sample = matfile('../../data/paris/features_sample_100k.mat');
Xval = data_sample.val_coefs;
Yval = data_sample.val_f;
Nval = data_sample.val_norm;

% properties
k = 5;

% data normalize
X_mean = mean(Xtr, 2);
X_std = std(Xtr, [], 2);
std_nz = find(X_std ~= 0);
%X_std(X_std == 0) = 1;
Xtr = Xtr(std_nz, :);
Xval = Xval(std_nz, :);
X_mean = X_mean(std_nz);
X_std = X_std(std_nz);
Xtr = (Xtr - repmat(X_mean, 1, size(Xtr, 2)))./repmat(X_std, 1, size(Xtr, 2));
Xval = (Xval - repmat(X_mean, 1, size(Xval, 2)))./repmat(X_std, 1, size(Xval, 2));
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

val_error = zeros(1, size(Xval, 2));
val_error_k = zeros(k, size(Xval, 2));
for i = 1:size(Xval, 2)
%for i = randperm(size(Xval, 2))
%for i = [51109, 66610, 70565, 52147, 65528]
  %repS = adprintf({}, [num2str(i), '/', num2str(size(Xval, 2))]);
  YvalDN = Yval(:, i).*Nval(:, i);
  dist = Xtr - repmat(Xval(:, i), 1, size(Xtr, 2));
  %dist = sum(dist.^2, 1);
  dist = sum(dist.*(Cinv*dist), 1); %#ok<MINV>
  Yoval = zeros(2, k);
  for kk = 1:k
    [minDist, idx] = min(dist);
    minDist = sqrt(minDist);
    dist(idx) = Inf;
    Yoval(:, kk) = Ytr(:, idx).*Ntr(:, idx);
    val_error_k(kk, i) = norm(YvalDN - Yoval(:, kk));
    fprintf([num2str(i), '/', num2str(kk), ': %5.2s %3.0f\n'], minDist, val_error_k(kk, i));
  end
  Yoval = mean(Yoval, 2);
  val_error(i) = norm(YvalDN - Yoval);
  fprintf([num2str(i), ': %3.0f\n'], val_error(i));
  %repS = rmprintf(repS);
end

[pathstr, name, ~] = fileparts(mfilename('fullpath'));
save([pathstr, '/results/', name, '-' datestr(now, 'yymmdd-HHMMSS') '.mat'], 'k', 'val_error', 'val_error_k');
clear X_mean X_std Xtr Xval Ytr Yval YvalDN Ntr Nval val_error val_error_k repS dist k kk i idx Yoval name pathstr;