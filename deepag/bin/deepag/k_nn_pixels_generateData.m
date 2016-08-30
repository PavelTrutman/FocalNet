% Script to generate syntetic data to k-NN testing in pixel space
% 
% Pavel Trutman
% INRIA, 2016

valSize = 500;
overMarginCoef = 1.05;

v = [repmat([-6/7; -4/7; -2/7; 0; 2/7; 4/7; 6/7], 2, 1); -repmat([-6/7; -4/7; -2/7; 0; 2/7; 4/7; 6/7], 2, 1)];

corr_tr.u = v;
corr_val.u = zeros(28, valSize);
for i = 1:valSize
  corr_val.u(:, i) = v + unifrnd(-1/7*overMarginCoef, 1/7*overMarginCoef, 28, 1);
  % permute the correspondences
  shuffle = randperm(7)';
  corr_val.u(:, i) = corr_val.u([shuffle + 0; shuffle + 7; shuffle + 14; shuffle + 21], i);
end

save('../../data/paris/correspondences_synteticKNN.mat', 'corr_tr', 'corr_val', '-v7.3');