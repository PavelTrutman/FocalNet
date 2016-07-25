% im=digicirc(r[,mr]) - Binary image of a digital circle
%
% r  ... radius
% mr ... window radius
% im ... image of the circle

% T. Pajdla, pajdla@neovision.cz, 2005 Aug 4
function im=digicirc(r,mr)

if nargin<2
    mr=r;
end

D  = 2*mr+1;
x  = -mr:mr;
A  = ones(D,1) * (x.^2) + (x.^2)' * ones(1,D);
im = ~ (A>=r^2);

return
