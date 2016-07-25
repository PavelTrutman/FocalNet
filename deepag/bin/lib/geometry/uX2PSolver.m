% P = uX2PSolver(X[,P]) - Absolute Camera Pose Solver
%
% P        = Camera P matrix celarray
% X(1:3,:) = 3 x 3 points in image 1
% X(4:6,:) = 3 x 3 3D points

% T. Pajdla, pajdla@cmp.felk.cvut.cz, 2015-09-09
function P = uX2PSolver(X,P)
if nargin>0
    if nargin<2 % solve
        P = {uX2P(X(1:3,:),a2h(X(4:6,:)))};
    end
else % unit tests
    P = true;
end
