% Script to generate syntetic data for fundamental 
% 
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016
function fund_generateData()
% generate 2 matrices: calibration, rotation and translation
% generate valSize random 8-tuples of 3d points.
% save valSize 8-tuples of correspondences 
valSize = 5;
overMarginCoef = 1.05;

v2K(unifrnd([5 0 0], [15 10 10]))

v=zeros(8,3);
corr.u = zeros(28, valSize);
for i = 1:valSize
  points = v + unifrnd(-1000, 1000, 8, 3)
  % permute the correspondences
  shuffle = randperm(8)';
  for j = 1:3
      points(:,j)=points(shuffle,j);
  end
  %corr.u(:, i) = 
end

% save('../../data/paris/correspondences_syntetic.mat', 'corr_tr', 'corr_val', '-v7.3');
end

function K = v2K(v)
% v is in the form {focal_length p1 p2}
K=[v(1) 0 v(2); 0 v(1) v(3); 0 0 1];
end
function K = X_(u)
K=[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]
end