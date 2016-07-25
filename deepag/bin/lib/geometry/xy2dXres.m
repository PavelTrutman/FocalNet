% [e,J] = xy2dXres(x,y,H) - Target points homography residuals 
%
% x     ... 2 x 1 2D target point 
% y{i}  ... i-th view 3 x n 2D points
% H{i}  ... i-th view homography
% e     ... residuals

% (c) T.Pajdla, www.neovision.cz, Feb 2 2006
function e = xy2dXres(x,y,H)
for i=1:length(H)
    % 'y-y' error: e{i} = xy2Hres(H{i},x,y(:,i),'y-y');
    xx(:,i) = h2a(H{i}*a2h(x));
    yy(:,i) = h2a(y(:,i));
    e(i) = sqrt(sum(yy(:,i)-xx(:,i)).^2)';
end
return
