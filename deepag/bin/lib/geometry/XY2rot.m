% R = XY2rot(X1,X2) - optimal rotation registration
%
% R = arg min ||X2 - r*X1||
%
% X1 ... 3 x n points
% X2 ... 3 x n points
% R  ... 3 x 3 rotation

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-08-09
function R = XY2rot(X1,X2)
% direction correlation matrix
K = zeros(3,3);				
for i=1:size(X1,2)						
  K = K + X2(:,i)*X1(:,i)';	
end			
[U,~,V] = svd(K);
S = diag([1 1 det(U)*det(V)]);
R = U*S*V';
