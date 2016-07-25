% [l,r] = l2Dfit(y) - Line fit
% y ... 2 x n points
% l ... 3 x 1 line coordinates
% r ... 1 x n residuals
% z ... 4 x n [the closest points z on the line; y]
% See l2Dres

% (c) T.Pajdla, www.neovision.cz, Sep 14, 2005
function  [l,r,z] = l2Dfit(y,~)
    my = mean(y,2);
    z  = y - my*ones(1,size(y,2));
    C  = z*z';
    [v,e] = eig(C);
    l  = [v(:,1)' -v(:,1)'*my]';
    % l  = l/sqrt(sum(l(1:2).^2)); not needed, v's are normalized
    if nargout>1
        r = l2Dres(l,y);
    end
    if nargout>2
        z = [y-l(1:2)*(l'*a2h(y));y];
    end
return