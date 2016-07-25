% [R,x,e] = ux2rd(u,X,R,rc,opts) - radial distortion estimation from homographies
% 
% u ... u{i} = 2 x n projections if the i-th target
% X ... 2 x n coordinates of the target
% R ... radial distortion descriptor
%       R.x0  ... 2 x 1 distortion center ([0;0] implicit if field r.x0 missing)
%       R.p   ... 1 x k ... k radial distortion coefficients
%       R.m   ... model type
%                 1 ... polynomial (implicit if field r.m missing)
%                       c = 1 + p(1)*r.^2 + p(2)*r.^4 + ...
%                 2 ... the division model
%                       r = r' * 1/(1+a*r'^2) reverse by r' = (1-sqrt(1-4*a*r^2))/(2*a*r)
% rc ... rectangle ~ [uMin uMax vMin vMax]
% e ... residuals
% x ... x{i} = undistorted points 
%
% See also uL2rdres, ur2u 

% (c) T.Pajdla, cmp.felk.cvut.cz, September 15, 2007
function [R,x,e,z] = ux2rd(u,X,R,rc,opts)
if nargin<5
    opts  = optimset('Display','off','TolX',1e-5,'MaxIter',100);
end
if nargin<4
    rc = [];
end
if ~isfield(R,'A')
    if isempty(rc)
        R.A = eye(2);
    else
        R.A = 4*sqrt(2)/((rc(2)-rc(1))^2+(rc(4)-rc(3))^2)*eye(2);
    end
end
% p = [R.x0(1) R.x0(2) sign(R.p).*abs(R.p).*10.^(9*(1:length(R.p)))];
p = [R.x0(1);R.x0(2);R.p];
m = R.m;
if opts.MaxIter>0
    p = lsqnonlin(@ux2rdres,p,[],[],opts,R,u,X); % optimize
    R.x0 = p(1:2); R.x0 = R.x0(:);
    % R.p = sign(p(3:end)).*abs(p(3:end))./10.^(9*(1:length(p)-2));
    R.p = p(3:end);
end
if nargout>1 % corrected points
    for i=1:length(u)
        x{i} = ur2u(u{i},R); % undistorted points
    end
end
if nargout>3 % residuals
    for i=1:length(u)
        [e{i},z{i}] = ux2rdres(p,R,u(i),X);
    end
else
    for i=1:length(u)
        e{i} = ux2rdres(p,R,u(i),X);
    end    
end
