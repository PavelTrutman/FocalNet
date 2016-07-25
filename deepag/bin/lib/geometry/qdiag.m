% [D,H] = qdiag(Q[,s]) - quadric diagonalization
% 
% Q ... 4x4 quadric
% s ... 0 or missing = diagonalization by an orthonormal H
%       1            = diagonalization by a Euclidean transform
%                      H = [R T; 0 0 0 1], R orthonormal
%                      R is chose such that x,y,z coordinates are assined
%                      to the quadric axes with their increasing lengths             

% H ... 4x4 change of the basis D = H'*Q*H

% (c) T.Pajdla, www.neovision.cz, Oct 28, 2004
function [D,H] = qdiag(Q,s)

Q = (Q+Q')/2;
Q = Q/norm(Q);

if nargin<2
    s = 0;
end

if s
    [V,e]   = eig(Q(1:3,1:3));
    [ev,ei] = sort(diag(e));
    d       = diag(ev(end:-1:1));
    V       = V(:,ei(end:-1:1));
    V       = V/det(V);
    T       = V'*inv(Q(1:3,1:3))*Q(1:3,4);
    H       = [V -V*T; 0 0 0 1];
    D       = H'*Q*H;
else
    [H,D] = eig(Q);
end




