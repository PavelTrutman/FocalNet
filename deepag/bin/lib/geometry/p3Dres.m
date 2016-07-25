% r = p3Dres(l,y) - Plane residuals
% y ... 3 x n points
% l ... 4 x 1 plane coordinates
% r ... 1 x n residuals
% 
% See p3Dfit

% (c) T.Pajdla, pajdla@cvut.cz, 2015-01-03
function  r = p3Dres(p,y,~,~)
    r  = p'*[y ; ones(1,size(y,2))];
return