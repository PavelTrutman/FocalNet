% N = null3x3r2(A) - Null space of 3x3 matrix A of rank 2 (faster than null)
% 
% A = 3x3 real matrix
% N = null space, i.e. A*N = 0

% pajdla@cmp.felk.cvut.cz
function N = null3x3r2(A)

A = A/max(abs(A(:))); %faster than sqrt(sum(A(:).^2)); % normalize A
N = [+A(1,2)*A(2,3)-A(1,3)*A(2,2) +A(2,2)*A(3,3)-A(2,3)*A(3,2) +A(3,2)*A(1,3)-A(3,3)*A(1,2)
     -A(1,1)*A(2,3)+A(1,3)*A(2,1) -A(2,1)*A(3,3)+A(2,3)*A(3,1) -A(3,1)*A(1,3)+A(3,3)*A(1,1)
     +A(1,1)*A(2,2)-A(1,2)*A(2,1) +A(2,1)*A(3,2)-A(2,2)*A(3,1) +A(3,1)*A(1,2)-A(3,2)*A(1,1)];
n = sum(N.^2); 
[n,i] = max(n); % the largest column of N 
N = N(:,i)/sqrt(n); % normalize
