% A = E2gAM(E) - Edge list E to graph non-symmetric adjacency matrix A 
%
% E = n x 2(3) matrix [[v1,v2,w12], ...] with vertice v1, v2 and weight w12 
%     wij are set to one for n x 2 matrix

% 2015-09-06 pajdla@cmp.felk.cvut.cz
function A = e2gam(E)
if size(E,2)<3
    E(:,3) = ones(size(E,1),1);
end
n = max(max(E(:,1:2)));
% A = sparse(n,n);
% A(sub2ind(size(A),E(:,1),E(:,2))) = E(:,3); % ineefficient
A = sparse(E(:,1),E(:,2),E(:,3),n,n);
A = A-diag(diag(A));
