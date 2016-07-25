% e = ux2rdcres(c,u,x) - radial distortion center estimation residuals
% 
% u     ... u{i} = 2 x n projections of the i-th target 
% x     ... x{i} = 2 x n calibration taget points
% c     ... 2 x 1 distortion center
% e     ... distance of x mapped by the best homography from the lines
%           l = [c]_x u
%
% Assumes that x and u were normalized.
%
% See also ux2rdc

% (c) T.Pajdla, cmp.felk.cvut.cz, May 1, 2006
function e = ux2rdcres(c,u,x)

% error('not finished ...');

C = a2h(c); % the distortion center
C = C/norm(C,2);
xc = xx(C); % vector product skew-symmetric matrix
xc = xc{1};

for i=1:length(u)
    % Construct the linear system to get the homography projecting x to
    % lines l = [c]_x u, i.e. solve for H in 
    %            l'*H*X  = 0
    % [u v w]*H*[x y z]' = 0
    % u*(H11*x+H12*y+H13*z)+v*(H21*x+H22*y+H23*z)+w*(H31*x+H32*y+H33*z)= 0
    % [u*x u*y u*z v*x v*y v*z w*x w*y w*z]*[H11 H12 H13 H21 ... H33]' = 0 
    %                                                           A * h  = 0
    % A = [u*x u*y u*z v*x v*y v*z w*x w*y w*z]
    % A = l([1 1 1 2 2 2 3 3 3],:)'.*X([1 2 3 1 2 3 1 2 3],:)'
    l = xc*a2h(u{i}); % the lines of the pencil
    l = l./([1;1;1]*vnorm(l));
    X = a2h(x{i});
    X = X./([1;1;1]*vnorm(X));
    A = l([1 1 1 2 2 2 3 3 3],:)'.*X([1 2 3 1 2 3 1 2 3],:)';
    [U,D,V] = svd(A);
    % get the "fundamental matrix"
    H1 = reshape(V(:,end-3),3,3);
    H2 = reshape(V(:,end-2),3,3);
    H3 = reshape(V(:,end-1),3,3);
    H4 = reshape(V(:,end-0),3,3);
    % use the constraint on the distortion center
    M = [(H1*C) (H2*C) (H3*C) (H4*C)];
    [U,D,V] = svd(M);
    a = V(:,end);    
    H{i} = a(1)*H1'+a(2)*H2'+a(3)*H3'+a(4)*H4';
    m{i} = l'*H{i}; % transformed pencil of lines
    m{i} = m{i}./(vnorm((m{i}(:,1:2))')'*[1 1 1]); % unit line normal vector
    e{i} = sum(m{i}'.*a2h(x{i})); % residual vectors
end
%% reshape the residuals into a vector
e = cat(2,e{:});
%me=mean(e(:).^2);
e = e(:);
