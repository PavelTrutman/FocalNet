% e = qres(Q,X) - evaluate meaningful residuals for a quadric
% 
% X ... 3xN points
% Q ... 4x4 quadric
% e ... residuals

% (c) T.Pajdla, www.neovision.cz, Oct 25, 2004
function e = qres(Q,X)

X       = a2h(X);
P       = Q*X;
n       = [1;1;1;1]*sqrt(sum(P(1:3,:).*P(1:3,:)));
P       = P./n;
e       = sum(P.*X)/2;




