% r = l2Dres(l,y) - Line residuals
% y ... 2 x n points
% l ... 3 x 1 line coordinates
% r ... 1 x n residuals
% 
% See l2Dfit

% (c) T.Pajdla, www.neovision.cz, Sep 14, 2005
function  r = l2Dres(l,y,p1,p2)
    r  = l'*[y ; ones(1,size(y,2))];
return