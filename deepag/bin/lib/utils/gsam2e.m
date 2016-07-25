% E = gSAM2E(A) - Graph symmetric adjacency matrix A to edge list E

% 2015-05-10 pajdla@cmp.felk.cvut.cz
function E = gSAM2E(A)

[i,j] = find(triu(A));
E = [i j];