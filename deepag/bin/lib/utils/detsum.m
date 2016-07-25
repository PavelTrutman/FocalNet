%function [d,m] = detsum(a,A) - determinant of a matrix pencil
%
% a ... pencil variables a(i)   = 0 ... free variable
%                        a(i) =\= 0 ... a(i) evaluated
% A(:,:,i) ... matrices of the pencil 
% d ... coefficients of the polynomial det(sum_i a(i)*A(:,:,i))
%       with monomials of free variables
% m ... multidegrees of monomials
%

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 25 Mar 2012
function [D,M] = detsum(a,A)

n = length(a); % length of the pencil
N = size(A,1); % size of NxN matrices
if false % This naive implementation works up to length(a)*size(A,1) < 10
    m = 1:n; % varible ids
    m = perms(repmat(m,1,N)); % possible variables
    m = m(:,1:N); % and their combinations up to deg N
    m = unique(m,'rows')';
else
    m = VChooseKRO(1:n,N)'; % k-tuples from 1:n with repetitionand ordering
end
v = a(m); v(isnan(v))=1; v = prod(v,1); % coefficients from bound variables
i = (n+1).^[(n-1):-1:0]; % numbering base 
i = i(m); % numbers
j = i.*isnan(a(m)); % only free variables are distinct
j = sum(j); % monomial index
[j,mi] = sort(j,'descend'); % sorted monomial index
m = m(:,mi); % sorted monomials
v = v(:,mi); % sorted coeffcients from bound variables
j = cumsum([1 diff(-j)>0]); % recode to get a linear index starting from 1
d = zeros(1,size(m,2)); % compute coefficients from determinants
D = zeros(1,length(unique(j))); % coefficients of the polynomial
g = diag(isnan(a)); % multidegree marix of free varibles
M = zeros(n,length(unique(j))); % monomials of the polynomial
for k=1:size(m,2)
    B = zeros(size(A,1),size(A,2)); % construct submatrix for monomial m(t,k)
    for t=1:N 
        B(:,t) = A(:,t,m(t,k));
    end
    d(k) = det(B); % determinant of the submatrix
    D(j(k)) = D(j(k)) + v(k)*d(k); % coeffcients of the final polynomial
    M(:,j(k)) = sum(g(:,m(:,k)),2); % monomilas of the final polynomial
end





         
         
         