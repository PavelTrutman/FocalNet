% y = E5ptResImPoints(x) - 5pt Minimal Relative Pose of Two Cameras - ransac residual evaluation data selection 
%
% x = struct 
%     x.x  = [x1;x2]   ... stacked homogeneous coordinates of image matches
%     x.iK = [iK1;iK2] ... stacked inverses of camera calibration matrices (jssumes triangular iKk)
% y = x.x

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-07-11
function y = E5ptResImPoints(x)
y = x.x;