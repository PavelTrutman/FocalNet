% R = xy2R3(x,y) - 3D rotation matrix R s.t. y = R x
%
% x,y ... column vectors
% R   ... rotation matrix s.t. y = R x

% (c) T. Pajdla, pajdla@neovision.cz
% 1 May 2008
function R = xy2R3(x,y)

x = x*diag(1./vnorm(x));
y = y*diag(1./vnorm(y));
switch size(x,2)
    case 1
        v = cross(x,y);
        ct = x'*y;
        st = norm(v);
        v = v / st;
        vt = 1-ct;
        R = [ ct	  -v(3)*st  v(2)*st
              v(3)*st  ct      -v(1)*st
             -v(2)*st  v(1)*st  ct	];
        R = v*v'*vt+R;
    case 2
        x1 = x(:,1);
        x2 = xx(x1)*x(:,2); x2 = x2/norm(x2);
        x3 = xx(x2)*x1;
        X = [x1 x2 x3];
        y1 = y(:,1);
        y2 = xx(y1)*y(:,2); y2 = y2/norm(y2);
        y3 = xx(y2)*y1;
        Y = [y1 y2 y3];
        R = Y*X';
    otherwise
        error('wrong size of input')
end

