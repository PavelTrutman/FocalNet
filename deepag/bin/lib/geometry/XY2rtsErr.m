% e = XY2rtsErr(E,X)
%
% e = ||X{2} - (s*r*X{1} + t)||
%
% X1 ... 3 x n points
% X2 ... 3 x n points
% E.r  ... 3 x 3 rotation
% E.t  ... 3 x 1 translation
% E.s  ... 1 x 1 scale
% e  ... errors

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 6 Jan 2009
function e = XY2rtsErr(E,X,cnstr)

e = sqrt(sum( (X(4:6,:)-(E.s*E.r*X(1:3,:)+E.t*ones(1,size(X,2)))).^2 ));
