% [P,in,e] = uX2PRFit(u,X,C,op) - Absolute camera pose from one image
%
% u = 2 x n image projections
% C = camera calibration structure, see X2u.m
% op = parameters
%      RANSAC (see, e.g. ransacfit with P4ptP4Pf model):
%      op.maxres = 2*sN;  % maximal residual in pixels = 2*max image noise
%      op.smpls  = 0;     % minimal sample size
%      op.smpln  = 1000;  % number of samples
%      op.doLO   = false; % do Local Optimization
%      op.cnstr  = [];    % constraints
% P = 3 x 4 camera projection matrix
% in = inlier selector
% e = reprojections error
% in = inlier selector 

% T. Pajdla, pajdla@cmp.felk.cvut.cz, 2015-09-08
function [P,in,e,x,f] = uX2PRFit(u,X,C,op)
if nargin>0
    xiK.x  = [a2h(u); X];% stacked uncalibrated image points & 3D points
    xiK.K = C.K; % camera calibration matrix
    xiK.iK = xiK.K\eye(3); % inverse of the camera calibration matrix
    [P,in] = ransacfit(xiK,op.solver,op.maxres,op.smpls,op.smpln,op.doLO,op.cnstr); % run ransac
    if ~isempty(in)
        e = PerspRepErr(P,xiK.x,xiK.K); % reprojection errors
        P = xiK.K*P; % get P for the original image coordinates
    else
        e = [];
    end
    f = []; % no figures will be plot
else % unit tests
    K = [1000 0 500;0 1000 500;0 0 1];
    R = a2r([1;1;1],pi/20);
    C = [0;0;-5];
    P = KRC2P(K,R,C);
    X = rand(3,1000)-0.5;
    u = h2a(X2u(X,P));
    op = struct('maxres',1,'smpls',0,'smpln',100,'doLO',false,'cnstr',[]);
    [Pr,in,e] = uX2PRFit(u,X,struct('K',K),op);
    test = acos(abs((P(:)/vnorm(P(:)))'*(Pr(:)/vnorm(Pr(:)))));
    P = test < [1e-7];
end

