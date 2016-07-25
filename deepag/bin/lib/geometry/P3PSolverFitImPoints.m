% y = P3PSolverFitImPoints(x) - 3pt Minimal Absolute Camera Pose - ransac calibrated fit data computation
%
% x = struct 
%     x.x  = homogeneous coordinates of image matches
%     x.iK = inverse of camera calibration matrix
% y = iK*x.x

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-07-11
function y = P3PSolverFitImPoints(x)
y = [x.iK*x.x(1:3,:)/x.iK(3,3); x.x(4:end,:)];