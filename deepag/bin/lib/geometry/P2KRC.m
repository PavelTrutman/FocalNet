% [K,R,C] = P2KRC(P) - Camera matrix P decomposition to K, R, C: P = K*R*[I -C]
%	 
%	P	    = camera calibration matrix ; size = 3x4
%	K	    = matrix of internal camera parameters
%	R	    = rotation matrix
%	C	    = camera center

% Tomas Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-08-26
function [K,R,C] = P2KRC(P)

if all(isfinite(P(:)))
    d = det(P(1:3,1:3));
    if abs(d)>100*eps
        P = P*sign(d);
    end
    B     = P(1:3,1:3);
    C     = -pinv(B)*P(:,4);
    [K,R] = kr3(B);
else % to work for nan matrices
    K = triu(nan(3)); R = nan(3); C = nan(3,1);
end
return  
