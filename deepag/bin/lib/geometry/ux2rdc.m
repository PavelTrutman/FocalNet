% c = ux2rdc(u,x[,mth,opts]) - radial distortion center
% 
% u     ... u{i} = 2 x n projections if the i-th target 
% x     ... 2 x n calibration taget points
% mth   ... 'svd' = svd fit (implicit)
%       ... 'bnd' = svd + bundle adjustment
% opts  ... lsqnonlin options, see optimset
% c     ... 2 x 1 distortion center
%
% Method by R.Hartley, S.B.Kang, ICCV 2005.
%
% See also ux2rdcres, ur2u 

% (c) T.Pajdla, cmp.felk.cvut.cz, May 1, 2006
function [c,F] = ux2rdc(u,x,mth,opts)
if nargin<4 && nargin>2
    opts  = optimset('Display','none','TolX',1e-5,'MaxIter',100);
end
if nargin<3
    mth = 'svd';
end
%% normalize points
Ax = x2nx(cat(2,x{:})); % normalization matrix for x
Au = x2nx(cat(2,u{:})); % normalization matrix for u
if 0
    Ax = eye(3);
    Au = eye(3);
end
for i=1:length(u) % normalize
    xn{i} = h2a(Ax*a2h(x{i})); 
    un{i} = h2a(Au*a2h(u{i})); 
end
for i=1:length(un)
    %% construct the linear system to get the "fundamental matrix"
    % A = [xu xv xw yu yv yw zu zv zw]
    uu = a2h(un{i});
    uu = uu([1 2 3 1 2 3 1 2 3],:)';
    xx = a2h(xn{i});
    xx = xx([1 1 1 2 2 2 3 3 3],:)';
    A = xx.*uu;
    if 0
        subfig(3,4,9);
        plot(vnorm(A'));
        title('vnorm(A'')');
    end
    if 0
        A = A./(vnorm(A')'*ones(1,size(A,2))); % balance A
    end
    % use eig instead of [U,D,V] = svd(A) to get it faster
    [V,D] = eig(A'*A);
    % get the "fundamental matrix"
    F{i} = reshape(V(:,1),3,3)';
end
%% get the distortion center
% F's are balanced already
f = cat(1,F{:}); 
% use eig instead of [U,D,V] = svd(f) to get it faster
[V,D] = eig(f'*f);
c = h2a(V(:,1));
%% bundle adjustment
if strcmp(mth,'bnd')
    % error('not finished ...');
    c = lsqnonlin(@ux2rdcres,c,[],[],opts,un,xn); % optimize
end
% undo the normalization
c = h2a(inv(Au)*a2h(c));
