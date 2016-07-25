% [q,x,V] = QP2q(Q,P[,N]) - Intersection of a quadric Q with a plane P
% 
% Q ... 4x4 quadric
% P ... 4x1 plane
% q ... 3x3 conic in the plane P
% V ... change of coordinates to make P into x=0
% x ... 4xN points on the intersection conic (N points 100 pts implicit),
%       homogeneous coordinates

% (c) T.Pajdla, www.neovision.cz, Oct 25, 2004
function [q,x,V]= QP2q(Q,P,N)

if nargin<3
    N = 100;
end

Q = (Q+Q')/2;
P = P(:);
if norm(P(1:3))>eps % Euclidean coordinate system
    P  = P/norm(P(1:3));
    n3 = P(1:3);
    d  = P(4);
    if norm(cross(n3,[0;0;1]))>eps % plane normal in a general position
        n2 = cross(n3,[0;0;1]);
        n1 = cross(n2,n3);
        V  = [[n1';n2';n3'] [0;0;d]
               0 0 0               1];
    else     % plane normal already along z axis
        V  = [eye(3) [0;0;d]
               0 0 0       1];
    end
    P  = inv(V)'* P;         % make plane z=0
    q  = inv(V)'*Q*inv(V);   % change the coordinate system for Q
    qi = logical([1;1;0;1]); % choose the conic in z=0
    q  = q(qi,qi);
    x  = sampleq(q,N);
    if ~isempty(x)
        y  = [x(1:2,:)
              zeros(1,size(x,2))
              x(3,:)];
        x  = inv(V)*y;
    end
else % projective coordinate system
    P       = P/norm(P); 
    PP      = (P*P')/2;
    [V,e]   = eig(PP);
    [mv,mi] = max(e(:));
    mi      = mod((mi-1),4)+1;
    qi      = ones(4,1);
    qi(mi)  = 0;
    qi      = logical(qi);
    q       = V'*Q*V;
    q       = q(qi,qi);
end

