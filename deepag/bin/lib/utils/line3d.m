%h = line3d(X,L,pars) - single matrix 3 x N line plot in 3D
% 
% X    ... single matrix 3 x N
% L    ... 2 x M indices of lines [[start1; end1], ...] into X
%
% See also LINE

% pajdla@cmp.felk.cvut.cz, 2015-09-26
function h = line3d(varargin)

X = varargin{1};
L = varargin{2};
if size(X,1)==3
    h0 = line([X(1,L(1,:));X(1,L(2,:))],[X(2,L(1,:));X(2,L(2,:))],[X(3,L(1,:));X(3,L(2,:))],varargin{3:end});
elseif size(X,1)==2
    h0 = line([X(1,L(1,:));X(1,L(2,:))],[X(2,L(1,:));X(2,L(2,:))],varargin{3:end});
else
    h0 = line(varargin{1:end});
end
if nargout>0
    h=h0;
end





