% X = uP2X(u,C,P) - 3D points X on a calibration planes P seen by pixel u
% 
% P ... 4 x 1 plane vector
% C ... camera 
%       C.K  ... internal calibration matrix
%       C.E  ... camera euclidean coordinate system transformation
%       C.rp ... [a] parameter of the radial distortion division
%                 model r = r' * 1/(1+a*r'^2)
% u ... 3 x n image points, affine coords for finite points
% X ... m x n x 4 homogeneous coordinates of 3D point X

% (c) T.Pajdla, www.neovision.cz, Oct 5, 2004
function X = uP2X(u,C,P)

L = u2L(u,C);
% for i=1:size(L,3)
%     X(:,i) = squeeze(L(:,:,i))*P';
% end
% A faster implementation of the above
S(1,:,:) = (ones(size(L,3),1)*P)';
S(2,:,:) = S(1,:,:);
S(3,:,:) = S(1,:,:);
S(4,:,:) = S(1,:,:);
T        = L.*S;
X        = squeeze(sum(T,2));
w = abs(X(end,:))>eps;
X = h2a(X,w);


