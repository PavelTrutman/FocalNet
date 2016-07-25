% [h,X] = cgplot(A[,X,lc,lb]) - gplot with vertices X or uniformly dostributed on a unit circle
%
% A = n x n graph adjacency matrix or 2 x m edges of the graph
% X = 2(s) x n coordinates of vertices (on circle if missing)
% lc = {line color, line style}, {'b','-'} if missing
% lb = vertex labels
% h = handles returned by line

% 2015-05-10 pajdla@cmp.felk.cvut.cz
function [h,X] = cgplot(A,X,lc,lb)

if nargin<2 || isempty(X)
    n = max(size(A));
    t = (0:2*pi/n:2*pi*(n-1)/n);
    X = [cos(t); sin(t)];
end
if nargin<3
    lc = {'b','-'};
end
if isempty(lc)
    lc = {'b','-'};
end
if size(A,1)~=2
    e = gam2e(A)';
else
    e = A;
end
h = line3d(X,e,'Color',lc{1},'LineStyle',lc{2});
if nargin>3
    text3d(1.05*X,lb);
end