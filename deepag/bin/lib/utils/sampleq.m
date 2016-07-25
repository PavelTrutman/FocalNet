% x = sampleq(q[,N]) - sample conic q
% 
% q  ... 3x3 conic
% N  ... number of points on the conic (100 inplicit)
% x  ... 3xN points on the conic (homogneous coords, 
%        use h2a to get affine coords)
%        = [] if no real points

% (c) T.Pajdla, www.neovision.cz, Oct 25, 2004
function [x,e,V] = sampleq(q,N)

if nargin < 2
    N = 100;
end

q     = (q+q')/2;
[V,e] = eig(q);
if e(1,1)*e(2,2)<0
    U = [0 0 1;
         0 1 0;
         1 0 0];
else
    U = eye(3);
end
V  = V*U;
e  = U'*e*U;
e  = e/sign(e(1,1))/abs(e(3,3));
if e(3,3)>0
    x = [];
    return
end
fi = (0:N-1)*2*pi/N;
x  = [cos(fi)/sqrt(e(1,1))
      sin(fi)/sqrt(e(2,2))];
x  = a2h(x);
x  = V*x;
