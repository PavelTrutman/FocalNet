% e = EGeomErr(E,X,iK) - 5pt Minimal Relative Pose of Two Cameras by D.Nister - Error
%
% E        = Fundamenral matrix(ces) as E(:,:), E(:,:,i) or E{i}
% X(1:2,:) = 3 x n uncalibrated points in image 1 (= x1) 
% X(3:4,:) = 3 x n uncalibrated points in image 2 (= x2)
% iK       = [iK1;iK2] ... stacked inverses of camera calibration matrices (implicit identities)
% e        = error as e(i,:) = max(||iK2*x2,E*iK1*x1||,||iK1*x1,(iK2*l2)'*E||)

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-07-11
function e = EGeomErr(E,X,iK)
if ~isempty(E)
    % get camera calibration matrices
    if nargin<3
        iK1  = [1 0 0;0 1 0;0 0 1];
        iK2t = iK1;
    else
        iK1  = iK(1:3,:);
        iK2t = iK(4:6,:)';
    end
    if iscell(E)
        E = cat(3,E{:});
    end
    x1 = X(1:3,:);
    x2 = X(4:6,:);
    e = zeros(size(E,3),size(X,2)); % initialize errors
    for i=1:size(E,3) % for all models
        F = iK2t*E(:,:,i)*iK1; % get fundamental matrix
        l2 = F*x1; % epipolar line
        % l2 = l2./([1;1;1]*(sqrt(sum(l2(1:2,:).^2)))); % unit normal vector
        e2 = abs(sum(l2.*x2)./sqrt(sum(l2(1:2,:).^2))); % point-line distance 
        l1 = F'*x2; % epipolar line
        % l1 = l1./([1;1;1]*(sqrt(sum(l1(1:2,:).^2)))); % unit normal vector
        e1 = abs(sum(l1.*x1)./sqrt(sum(l1(1:2,:).^2))); % point-line distance
        e(i,:)  = max(e1,e2); % max of them 
    end
else
    e = [];
end