% [Q,e] = Qfit(X) - LSQ Quadric fit 
% 
% X ... 3xN points
% Q ... 4x4 quadric
% e ... residuals

% (c) T.Pajdla, www.neovision.cz, Oct 25, 2004
function [Q,e] = Qfit(X)

% B = [eye(3)/s -X0;0 0 0 1]
% Y = B*X

X0 = mean(X')';
Y  = X-X0*ones(1,size(X,2));
s  = std(sqrt(sum(Y.*Y)));
Y  = Y/s;
B  = [eye(3)/s -X0/s;0 0 0 1];
% 0 = Y'*O*Y
% 0 = (B*X)'*O*(B*X) = X'*B'*O*B*X => Q = B'*O*B

Y       = a2h(Y);
A       = [([1;1;1;1]*Y(1,:)).*Y
             ([1;1;1]*Y(2,:)).*Y(2:4,:)
               ([1;1]*Y(3,:)).*Y(3:4,:)
                     (Y(4,:)).*Y(4:4,:)];         
A       = diag([1 2 2 2 1 2 2 1 2 1])*A;
A       = A';
[u,d,v] = svd(A,0);
q       = v(:,end);
Q       = [q(1) q(2) q(3) q(4)
           q(2) q(5) q(6) q(7)
           q(3) q(6) q(8) q(9)
           q(4) q(7) q(9) q(10)];
Q       = B'*Q*B;
Q       = Q/norm(Q);
e       = qres(Q,X);

