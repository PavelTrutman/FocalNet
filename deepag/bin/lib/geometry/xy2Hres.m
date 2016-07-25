% e = xy2Hres(H,x,y,obj) - Homography residuals 
%
% H     ... 3 x 3 (9 x 1) homography matrix
% x     ... 3 x n 2D points 
% y     ... 3 x n 2D points
% obj   ... 'x-x' ... error in x
%                     e  = dy = sqrt(sum((h2a(inv(H)*y)-h2a(x)).^2));
%           'y-y' ... error in y
%                     e  = dx = sqrt(sum((h2a(y)-h2a(H*x)).^2)); 
%           'x-y' ... mean odf the error in x + error y (implicit)
%                     e  = (dx + dy)/2;
%           'xmy' ... min of the error in x and error in y 
%                     e  = min(dx,dy);

% (c) T.Pajdla, www.neovision.cz, Nov 2 2005
function e = xy2Hres(H,x,y,obj)
    if any(size(H)-[3 3])
        H = reshape(H,3,3);
    end
    switch obj
        case 'x-x'
            e = sqrt(sum((h2a(inv(H)*y)-h2a(x)).^2))';
        case 'y-y'
            e = sqrt(sum((h2a(y)-h2a(H*x)).^2))';
        case 'x-y'
            e = (sqrt(sum((h2a(y)-h2a(H*x)).^2))+sqrt(sum((h2a(inv(H)*y)-h2a(x)).^2)))'/2;
        case 'xmy'
            e = min(sqrt(sum((h2a(y)-h2a(H*x)).^2)),sqrt(sum((h2a(inv(H)*y)-h2a(x)).^2)))';
        case 'xXy'
            e = max(sqrt(sum((h2a(y)-h2a(H*x)).^2)),sqrt(sum((h2a(inv(H)*y)-h2a(x)).^2)))';            
    end
return
