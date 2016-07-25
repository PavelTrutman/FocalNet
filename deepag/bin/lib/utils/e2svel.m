% e = e2svel(E) - Edge list E to symmetric vertex edge list
%
% E = n x 2 matrix [[v1,v2], ...] with vertice v1, v2 
% e = {[e1 e2 e3 ...], ...} cellaray of lists of edges ids incident to a vertex

% 2015-09-06 pajdla@cmp.felk.cvut.cz
function e = e2svel(E)

mV = max(E(:)); % maximal vertex id
e = cell(mV,1); % prepare a cellarray
for i=1:size(E,1)
    e{E(i,1)}(end+1) = i;
    e{E(i,2)}(end+1) = i;
end
