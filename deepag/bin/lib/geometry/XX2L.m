% L = XX2L(X,Y) - Plucker line matrix from 2 points
% 
% X ... 4 x 1 [x n] point[s]
% Y ... 4 x 1 point
% L ... 4 x 4 [x n] ray Plucker matrix[ces]

% (c) T.Pajdla, www.neovision.cz, Feb 3, 2005
function L= XX2L(X,Y)

if size(X,2)==1
    L = X*Y'-Y*X';
else
    L(:,1,:) = X * Y(1);
    L(:,2,:) = X * Y(2);
    L(:,3,:) = X * Y(3);
    L(:,4,:) = X * Y(4);    
    L(1,:,:) = squeeze(L(1,:,:)) - X * Y(1);
    L(2,:,:) = squeeze(L(2,:,:)) - X * Y(2);    
    L(3,:,:) = squeeze(L(3,:,:)) - X * Y(3);    
    L(4,:,:) = squeeze(L(4,:,:)) - X * Y(4);       
end
