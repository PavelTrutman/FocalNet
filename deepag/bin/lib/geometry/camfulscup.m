% [X,X1,X2] = camfulscup(P1,z1,P2,z2,[,tol,quick]) - perspective camera fulcra intesection
%
% P1, P2 = 3x4 perspective camea projection matrix P = K*R*[I|-C] with rotation R, pose C and
%
%          [1/bx s/f  u0/f] 
%     K  = [0    1/by v0/f]   
%          [0     0    1/f]   
%
%     determining the camera coordinate system \beta. 
% z1, z2 = [nz, fz]  near (nz) and far (fz) z clipping plane in space 
% X      = 3 x n vertices of the intersetion polyhedra (empty if empty)
% X1, X2 = 3 x 8 vertices of the fulcra 
% tol    = tollerance on intersection >-tol [1e-12]
% quick  = if true then exit asap there is asingle point inf the intersection discovered [false]

% Tomas Pajdla (pajdla@cmp.felk.cvut.cz)
% 2015-12-01
function [X,X1,X2] = camfulscup(P1,z1,P2,z2,tol,quick)
if nargin<6
    quick = false;
end
if nargin<5
    tol = 1e-12;
end
% fulcra
[X1,L1] = camfulcrum(P1,z1(1),z1(2));
[X2,L2] = camfulcrum(P2,z2(1),z2(2));
% intersect them
L = [L1;L2]; % all planes
ix = nchoosek(1:size(L,1),3); % all (220) possible triplets of different planes
empty = true; X = [];
for i=1:size(ix,1) % while no point found inside all planes
    n = L(ix(i,:),:); % 3 planes
    x = [h2a([det(n(:,[2 3 4]));-det(n(:,[1 3 4]));det(n(:,[1 2 4]));-det(n(:,[1 2 3]))]);1]; % intersection
    if all(L*x>-tol),
        empty = false;
        X(:,end+1) = x;
        if quick, break; end
    end
end
if ~empty
    X = X(1:3,:);
end
