% L = vlist(v) - Celarray of indices of occurences of different positive integer values in v
%
% v = a vector of positive integres (index)
% L = {[i1 i2 i3 ...], ...} cellaray of lists of occurenves of values in v L{i} = find(v==i) 

% 2015-09-06 pajdla@cmp.felk.cvut.cz
function L = vlist(v)

mV = max(v(:)); % maximal vertex id
L = cell(mV,1); % prepare a cellarray
for i=1:length(v)
    L{v(i)}(end+1) = i;
end
