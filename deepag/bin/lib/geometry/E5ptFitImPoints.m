% y = E5ptFitImPoints(x) - 5pt Minimal Relative Pose of Two Cameras - ransac calibrated fit data computation
%
% x = struct 
%     x.x  = [x1;x2]   ... stacked homogeneous coordinates of image matches
%     x.iK = [iK1;iK2] ... stacked inverses of camera calibration matrices (assumes triangular iKk)
% y = [iK1*x.x(1:3,:);iK1*x.x(1:3,:)]

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-07-11
function y = E5ptFitImPoints(x)

y = [x.iK(1:3,:)*x.x(1:3,:)/x.iK(3,3)
     x.iK(4:6,:)*x.x(4:6,:)/x.iK(6,3)];