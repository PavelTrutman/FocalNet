% n = vnorm(x) - columnwise vector 2 norm
% 
% n = sqrt(sum(x.^2))

% (c) T.Pajdla, cmp.felk.cvut.cz, May 2 2006
function n = vnorm(x)
    n = sqrt(sum(x.^2));
return
