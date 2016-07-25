% e = PerspRepErr(P,X,K) - Perspective reprojection error
%
% P        = Calibrated projection matrix [R -R*C], P, P(:,:,i) or P{i}
% K        = Camera calibration matrix
% X(1:3,:) = 3 x n image projections
% X(4:6(7),:) = 3(4) x n 3D points
% e        = image reprojection error

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-09-09
function e = PerspRepErr(P,X,K)
if ~isempty(P)
    if size(X,1)<7
        X = a2h(X);
    end
    if nargin<3
        K = eye(3);
    end
    if iscell(P)
        P = cat(3,P{:});
        if iscell(K)
            K = cat(2,K{:});
        else
            K = repmat(K,[1 1 size(P,3)]);
        end
    end
    e = zeros(size(P,3),size(X,2)); % initialize errors
    if size(P,3)<2
        e = vnorm(h2a(X(1:3,:))-h2a(K*P*X(4:7,:)));
    else
        for i=1:size(P,3)
            e(i,:) = vnorm(h2a(X(1:3,:))-h2a(K(:,:,i)*P(:,:,i)*X(4:7,:)));
        end
    end
else
    e = [];
end