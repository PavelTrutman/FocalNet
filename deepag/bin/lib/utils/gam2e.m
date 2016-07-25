% E = gAM2E(A) - Graph non-symmetric adjacency matrix A to edge list E

% 2015-05-10 pajdla@cmp.felk.cvut.cz
function E = gAM2E(A)

[i,j] = find(A);
E = [i j];