% [A,B] = xy2nxy(x,y[,df]) - Point normalization for multiple view geometry - 2 images
% 
% mean(A*x) = 0 & mean(B*y) = 0 & std(|A*x|) = sqrt(2) & std(|B*y|) = sqrt(2)
%
% x, y ... 2 x n points (if 3 x n, the last coord. is ignored)
% A, B ... 3 x 3 affinity matrix
% df   ... debug level flag, 0 or missing ... no debug info
%               
% R.Hartley, A.Zisserman. Multiple View Geometry. 2nd edition. Cambridge
% Press 2003.
%
% See also X2NX

% (c) T.Pajdla, www.neovision.cz, Nov 2 2005
function [A,B] = xy2nxy(x,y,df)

if nargin < 3
    df = 0;
end
if df > 0
    xo = x(1:2,:);
    yo = y(1:2,:);
end

x  = x(1:2,:);                  % ignore the last coords
x0 = mean(x')';                 % the centroid
x  = x-x0*ones(1,size(x,2));    % centered x
y  = y(1:2,:);                  % ignore the last coords
y0 = mean(y')';                 % the centroid
y  = y-y0*ones(1,size(y,2));    % centered y
k  = mean([sqrt(sum([x].^2)) sqrt(sum([y].^2))]);   % mean of lengths of x & y
A  = diag([sqrt(2)*[1 1]/k 1])*[eye(2) -x0; 0 0 1]; % construct the matrix
B  = diag([sqrt(2)*[1 1]/k 1])*[eye(2) -y0; 0 0 1]; % construct the matrix      

if df > 0
    subfig(3,4,1);
    % plot original plints
    plot(xo(1,:),xo(2,:),'.b','markersize',15);
    hold
    plot(yo(1,:),yo(2,:),'.r','markersize',10);    
    title('original points: b - x, r - y');
    subfig(3,4,2);
    xo = h2a(A*a2h(xo));
    yo = h2a(B*a2h(yo));    
    % plot transformed points
    plot(xo(1,:),xo(2,:),'.b','markersize',15);
    hold
    plot(yo(1,:),yo(2,:),'.r','markersize',10);    
    title('transformed points: b - x, r - y');
end
