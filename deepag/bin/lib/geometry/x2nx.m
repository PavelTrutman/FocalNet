% A = x2nx(x) - Point normalization for multiple view geometry
% 
% mean(A x) = 0  &  std(|A x|) = sqrt(2)
%
% x ... 2 x n points (if 3 x n, the last coord. is ignored)
% A ... 3 x 3 affinity matrix
%               
% R.Hartley, A.Zisserman. Multiple View Geometry. 2nd edition. Cambridge
% Press 2003.
%
% See also XY2NXY

% (c) T.Pajdla, www.neovision.cz, Nov 1 2005
function A = x2nx(x)

x  = x(1:2,:);                  % ignore the last coords
x0 = mean(x')';                 % the centroid
x  = x-x0*ones(1,size(x,2));    % centered x
d  = mean(sqrt(sum(x.^2)));     % mean of lengths of x
A  = diag([sqrt(2)*[1 1]/d 1])*[eye(2) -x0; 0 0 1]; % construct the matrix
      