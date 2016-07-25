% h = plotx2y(x,y[,s,xp,yp,vp]) - plot a difference fied from x to y
% 
% x,y ... 2 x n real matrix
% s   ... scale (implicit = 1)
% xp  ... x plot string ([] to ignore, implicit blue dot)
% yp  ... y plot string ([] to ignore, implicit green dot)
% vp  ... y-x plot string ([] to ignore, implicit red line)

% (c) T.Pajdla, pajdla@cmp.felk.cvut.cz
% 12 May 2010
function h = plotx2y(varargin)

if nargin<3, s = 1; else s = varargin{3}; end
if nargin<4, xp = []; else xp = varargin{4}; end
if nargin<5, yp = []; else yp = varargin{5}; end
if nargin<5, vp = []; else vp = varargin{6}; end

x = varargin{1};
y = varargin{2};
v = s*(y-x);

hd = ishold;
if isempty(xp) h0{1} = plot3d(x,'.b'); else h0{1} = plot3d(x,xp); end; hold on
if isempty(yp) h0{2} = plot3d(y,'.g'); else h0{2} = plot3d(y,yp); end
if size(x,1)<3
    lx = [x(1,:);x(1,:)+v(1,:)];
    ly = [x(2,:);x(2,:)+v(2,:)];    
    if isempty(vp) h0{3} = line(lx,ly,'color','r'); else h0{3} = line(lx,ly,'color',vp); end
else
    lx = [x(1,:);x(1,:)+v(1,:)];
    ly = [x(2,:);x(2,:)+v(2,:)];        
    ly = [x(3,:);x(3,:)+v(3,:)];    
    if isempty(vp) h0{3} = line(lx,ly,lz,'color','r'); else h0{3} = line(lx,ly,lz,'color',vp); end    
end
if ~hd, hold off; end

if nargout>0, h=h0; end





