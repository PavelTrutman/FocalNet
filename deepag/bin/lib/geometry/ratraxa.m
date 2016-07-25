% R = ratraxa(v,a) - a rational approximation of rotation aroung v by a
%
% v = rotation axis vector
% a = rotation angle in (-pi/2 pi/2)
% R = rotation matrix

% 2011-10-31, T. Pajdla, pajdla@gmail.com
function R = ratraxa(v,a)

tol = 10^6*eps;

v = v/norm(v);
[a(1),a(2)] = rat(tan(a),tol);
f = atan2(v(2),v(1));
t = atan2(norm(v(1:2)),v(3));
R = [a f t];
