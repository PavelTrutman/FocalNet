% d = XdX(x) - 2D point distance matrix
% 
% x ... 2 x n point coordinates
% d ... n x n point symmetric matrix of distances

% (c) T.Pajdla, www.neovision.cz, Sep 6, 2005
function d = XdX(x)

dm = ([repmat(x,[size(x,2),1])-repmat(x(:),[1,size(x,2)])]);
d  = sqrt(dm(1:2:end,:).^2+dm(2:2:end,:).^2);

