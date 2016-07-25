% d = XdL(X,L) - Distance of a point X from a line L
% 
% X ... 4 x n points
% L ... 4 x 4 line Plucker matrix
% d ... 1 x n distances

% (c) T.Pajdla, www.neovision.cz, Oct 5, 2004
function d = XdL(X,L)

[u,d,v] = svd(L);
PIn     = v(:,3:4);      
PIn     = [PIn./(ones(4,1)*sqrt(sum(PIn(1:3,:).^2)))]';
d       = sqrt(sum((PIn*X).^2));  


