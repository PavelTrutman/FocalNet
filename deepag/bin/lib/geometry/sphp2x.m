% x = sphs2x(p) - Spherical parametrization of a unit vector to the unit vector
% 
% sphs2x(p) = [cos(p(end))*sphs2x(p(1:end-1)) sin(p(end))]
% 
% p ... n-1 x 1 matrix [radians]
% x ... n x 1 matrix
%
% See also X2SPHS

% (c) T.Pajdla, www.neovision.cz, Nov 1 2005
function x = sphs2x(p)

x    = zeros(length(p),1);
x(1) = 1; 
for i=1:length(p)
    x(1:i+1) = [cos(p(i))*x(1:i); sin(p(i))];
end