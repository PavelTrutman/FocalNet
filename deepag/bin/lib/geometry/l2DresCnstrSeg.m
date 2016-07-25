% r = l2DresCnstrSeg(l,y,c) - Line residuals
% y ... 2 x n points
% l ... 3 x 1 line coordinates
% c ...[x1 x2], such that x1 resp. x2 is on the left resp. right 
%                                     from the line 
% r ... 1 x n residuals, Inf if constraints are not met
% 
% See l2Dfit

% (c) T.Pajdla, www.neovision.cz, Sep 14, 2005
function  r = l2DresCnstrSeg(l,y,c)
    r  = l'*[y ; ones(1,size(y,2))]; % residuals
    v  = (l'*[c;1 1]); % sidedness of points
    if v(1)*v(2)>0 % on different sides
        r = Inf*abs(r);
    end
return