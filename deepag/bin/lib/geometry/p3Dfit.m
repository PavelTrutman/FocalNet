% [l,r] = p3Dfit(y) - Plane fit
% y ... 3 x n points
% l ... 4 x 1 plane coordinates
% r ... 1 x n residuals
% z ... 6 x n [the closest points z on the line; y]
% See p3Dres

% (c) T.Pajdla, pajdla@cvut.cz, 2015-01-03
function  [p,r,z] = p3Dfit(y,~)
    my = mean(y,2);
    z  = y - my*ones(1,size(y,2));
    C  = z*z';
    [v,e] = eig(C);
    p  = [v(:,1)' -v(:,1)'*my]';
    % l  = l/sqrt(sum(l(1:2).^2)); not needed, v's are normalized
    if nargout>1
        r = p3Dres(p,y);
    end
    if nargout>2
        z = [y-p(1:3)*(p'*a2h(y));y];
    end
return