% P = P4ptP4Pf(X,P) - 4pt Minimal Ablsoute Pose Problem
%
% P        = 3x4 projection matrix
% X(1:3,:) = 3 x 4 3D points
% X(4:5,:) = 2 x 4 image projections

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-09-06
function P = P4ptP4Pf(X,P)

if nargin>0
    if all(size(X)==[5 4])
        [f,R,t] = P4Pf(X(4:5,:),X(1:3,:));
        P = cell(1,length(f));
        for i=1:length(f)
            P{i} = diag([f(i) f(i) 1])*[R(:,:,i) t(:,i)];
        end
    end
else % unit tests
    P = true;
end