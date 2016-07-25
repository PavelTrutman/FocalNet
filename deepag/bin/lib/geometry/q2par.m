% [p,A] = q2par(q) - Conic matrix q to conic parameters p = [x0 y0 a b th]
%
% q  = 3 x 3 conic matrix
% a  = major semiaxis
% b  = minor semiaxis
% th = rotation angle
% A = diagonalization matrix A'*q*A = diag([1/a^2 1/b^2 -1])

% T. Pajdla pajdla@ciirc.cvut.cz 2014-01-10
function [p,A] = q2par(q)

q = (q+q')/2; % symmetrize
q = q/max(abs(q([1 2 4 5]))); % normalize 2x2x submatrix for numerical conditioning
q = -q*sign(det(q)); % make the determinant negative
x  = - q(1:2,1:2) \ q(1:2,3); % the center
x0 = x(1); y0 = x(2);
T = [1 0 -x0;0 1 -y0;0 0 1]; % translation
Q = inv(T)'*q*inv(T); % move to the origin
Q = Q/abs(Q(3,3)); 
[v,e] = eig(Q(1:2,1:2)); % eigenvales/eigenvectors
s = v(1,2); c = v(1,1);  % the sine and cosine of the angle 
th = atan2(s,c); % the angle 
e = sign(e)*(1./sqrt(abs(diag(e))));
a = max(e); % major semi-axis
b = min(e); % minor semi-axis
p  = [x0;y0;a;b;th];
A = inv(T)*[c -s 0;s c 0;0 0 1];

