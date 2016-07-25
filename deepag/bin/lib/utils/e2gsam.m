% A = E2gSAM(E) - Edge list E to graph symmetric adjacency matrix A 

% 2015-09-06 pajdla@cmp.felk.cvut.cz
function A = e2gsam(E)

A = e2gam(E);
A = A+A'-diag(diag(A));
