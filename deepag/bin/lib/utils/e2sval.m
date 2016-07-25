% a = e2sval(E) - Edge list E to symmetric vertex adjacency list
%
% E = n x 2 matrix [[v1,v2], ...] with vertice v1, v2 
% a = {[v1 v2 v3 ...], ...} cellaray of lists of sjacent vertices 

% 2015-09-06 pajdla@cmp.felk.cvut.cz
function a = e2sval(E)

E = unique(E,'rows','sorted'); % remove multiplicities and sort rows to speedup the construction
mV = max(E(:)); % maximal vertex id
a = cell(mV,1); % prepare a cellarray
for i=1:size(E,1)
    a{E(i,1)}(end+1) = E(i,2);
    a{E(i,2)}(end+1) = E(i,1);
end
