% p = x2sphs(x) - Spherical parametrization of a unit vector
% 
% x2sphs(x) = [x2sphs(x(1:end-1))*acos(|x(1:end)|) asin(x(end))]
% 
% x ... n x 1 matrix
% p ... n-1 x 1 matrix
%
% See also SPHS2X

% (c) T.Pajdla, www.neovision.cz, Nov 1 2005
function p = x2sphs(x)

i = length(x)-1;
n = norm(x(1:i));
while n>eps
    p(i) = atan2(x(i+1),n);
    x    = x/n;
    i    = i - 1;    
    n    = norm(x(1:i));
end
