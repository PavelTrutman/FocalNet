% [X,L] = camfulcrum(P,z1,z2) - perspective camera fulcrum
%
% P = 3x4 perspective camea projection matrix P = K*R*[I|-C] with rotation R, pose C and
%
%          [1/bx s/f  u0/f] 
%     K  = [0    1/by v0/f]   
%          [0     0    1/f]   
%
%     determining the camera coordinate system \beta. 
% z1, z2 =  near and far clipping plane in space 
% X = 3 x 8 vertices
% L = 6 x 4 oriented planes (inside is positive)

% Tomas Pajdla (pajdla@cmp.felk.cvut.cz)
% 2015-12-01
function [x,p] = camfulcrum(P,z1,z2)
[K,R,C] = P2KRC(P); % camera parameters
K = K/K(3,3);
p_a = K(1:2,3); % principal point in \alpha
x = [0 2*p_a(1) 2*p_a(1) 0       
     0 0        2*p_a(2) 2*p_a(2)]; % outer image plane frame in \alpha
x = K\a2h(x); % outer frame in \epsilon
x = [z1*x z2*x]; % near and far clipping frames
x = [R' C]*a2h(x); % in \delta
p = [x2p(x(:,[1 2 3]))
     x2p(x(:,[6 5 8]))
     x2p(x(:,[6 2 1]))
     x2p(x(:,[2 6 7]))     
     x2p(x(:,[4 3 7]))     
     x2p(x(:,[5 1 4]))]; % planes 
 
function p = x2p(X)
  p = cr(X(:,2)-X(:,1),X(:,3)-X(:,2))';
  p = p/sqrt(sum(p.^2));
  p = [p -p*X(:,1)];
return

function z = cr(x,y)
  z = [x(2)*y(3)-x(3)*y(2);-x(1)*y(3)+x(3)*y(1);x(1)*y(2)-x(2)*y(1)];
return