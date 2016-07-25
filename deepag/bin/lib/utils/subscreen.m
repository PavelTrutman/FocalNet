% i = subscreen(r,c,rs,cs,is) - subfig helper: index i into r x c screen from index is into rs([1 2]) x cs([1 2]) subscreen range

% pajdla@cmp.felj.cvut.cz 2015-09-16
function i = subscreen(r,c,rs,cs,is)
if ~all(rs>0 & rs<=r) || ~all(cs>0 & cs<=c), error('rs, cs out of range'); end

sz = [diff(rs)+1 diff(cs)+1]; % size of the subscreen
is = mod(is-1,prod(sz)); % wrap index around & from 0 to prod(sz)-1
ris = floor(is/sz(2)); % row-1 in the subscreen
cis = mod(is,sz(2)); % column-1 in the subscreen
ri = rs(1)+ris; % row in the screen
ci = cs(1)+cis; % col in the screen
i = (ri-1)*c+ci; % index in the screen