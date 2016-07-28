% Script to generate syntetic data to k-NN testing in pixel space
% 
% Pavel Trutman
% INRIA, 2016

valSize = 500;
overMarginCoef = 1.05;

v = [repmat([-6/7; -4/7; -2/7; 0; 2/7; 4/7; 6/7], 2, 1); -repmat([-6/7; -4/7; -2/7; 0; 2/7; 4/7; 6/7], 2, 1)];

u_tr = v;
u_val = zeros(28, valSize);
for i = 1:valSize
  u_val(:, i) = v + unifrnd(-1/7*overMarginCoef, 1/7*overMarginCoef, 28, 1);
  % permute the correspondences
  shuffle = randperm(7)';
  u_val(:, i) = u_val([shuffle + 0; shuffle + 7; shuffle + 14; shuffle + 21], i);
end

save('pixels_syntetic.mat', 'u_tr', 'u_val', '-v7.3');