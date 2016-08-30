% A = x2nx(x[,mtd]) - Point normalization for multiple view geometry
%
% mtd = 'HZ'      ... mean(A x) = 0  &  std(|A x|) = sqrt(2) (implicit)
%     = '[-1,1]'  ... scaled to fit range(|A x|) = [-1 1]
%
% x ... 2 x n points (if 3 x n, the last coordinate is ignored)
% A ... 3 x 3 affinity matrix
%
% R.Hartley, A.Zisserman. Multiple View Geometry. 2nd edition. Cambridge
% Press 2003.
%
% See also X2NX

% (c) T.Pajdla, pajdla@cvut.cz, 2016-08-26
function A = x2nx(x,mtd)
if nargin>0
    if nargin<2
        mtd = 'HZ';
    end
    x  = x(1:2,:);                  % ignore the last coords
    switch mtd
        case 'HZ'
            x0 = mean(x')';                 % the centroid
            x  = x-x0*ones(1,size(x,2));    % centered x
            d  = mean(vnorm(x));     % mean of lengths of x
            A  = diag([sqrt(2)*[1 1]/d 1])*[eye(2) -x0; 0 0 1]; % construct the matrix
        case '[-1,1]'
            d  = vnorm(x); % lengths of x
            A = diag([(1/max(abs(x(:))))*[1 1] 1]);
    end
else
    xo = rand(2,100)+3*rand(1);
    subfig(3,4,1);
    % plot original points
    plot3d(xo,'.b'); title('original points');
    % HZ:  plot transformed points
    subfig(3,4,2);
    A=x2nx(xo,'HZ'); x = h2a(A*a2h(xo));
    plot3d(x,'.b'); title('HZ: transformed points');
    % scale: plot transformed points
    subfig(3,4,3);    
    A=x2nx(xo,'scale'); x = h2a(A*a2h(xo));
    plot3d(x,'.b'); title('scale: transformed points');    
    A = true;
end