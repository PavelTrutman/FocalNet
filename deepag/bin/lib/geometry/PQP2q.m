% [q,H] = PQP2q(Q,p,P) - Projection of the intersection of a quadric Q with a plane p by camera P
% 
% Q ... 4x4 quadric
% p ... 4x1 plane
% P ... 3x4 camera projection matrix
% q ... 3x3 conic in the plane P
% H ... change of coordinates to make p into w=0

% (c) T.Pajdla, www.neovision.cz, Oct 26, 2004
function [q,H] = QP2q(Q,p,P)

Q       = (Q+Q')/2;
p       = p(:)/norm(p);
[u,d,v] = svd(p');  
H       = inv(v)';
QQ      = inv(H)'*Q*inv(H);
PP      = P*inv(H);
A       = PP(:,2:4);
q       = inv(A)'*QQ(2:4,2:4)*inv(A);
