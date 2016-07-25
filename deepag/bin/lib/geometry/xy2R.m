% R = xy2R(x,y) - 2D rotation matrix R s.t. y = R x
%
% x,y ... column vectors
% R   ... rotation matrix s.t. y = R x

% (c) T. Pajdla, pajdla@neovision.cz
% 13 Aug 2007
function R = xy2R(x,y)

x = x*diag(1./vnorm(x));
y = y*diag(1./vnorm(y));
R = [ x'*y           [-y(2) y(1)]*x
    -[-y(2) y(1)]*x           x'*y];
