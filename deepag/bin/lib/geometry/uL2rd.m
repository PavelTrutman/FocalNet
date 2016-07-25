% [R,x,e,z] = uL2rd(u,L,R,rc,opts) - radial distortion estimation from lines
% 
% u ... u{i} = 2 x n projections if the i-th target
% L ... L{i}{j} = indices of points lying on a stright line in the tartget
% R ... radial distortion descriptor, see u2ru 
%       R.x0  ... 2 x 1 distortion center ([0;0] implicit if field r.x0 missing)
%       R.p   ... 1 x k ... k radial distortion coefficients
%       R.m   ... model type
%                 1 ... polynomial (implicit if field r.m missing)
%                       c = 1 + p(1)*r.^2 + p(2)*r.^4 + ...
%                 2 ... the division model
%                       r = r' * 1/(1+a*r'^2) reverse by r' = (1-sqrt(1-4*a*r^2))/(2*a*r)
%       R.T   ... normalization transform
% rc ... rectangle ~ image size
% x ... x{i} = undistorted points 
% e ... residuals
% z ... [the closet points on the lines; x]
%
% See also uL2rdres, ur2u 

% (c) T.Pajdla, cmp.felk.cvut.cz, July 13, 2007
function [R,x,e,z] = uL2rd(u,L,R,rc,opts)
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
% p = [R.x0(1) R.x0(2) sign(R.p).*abs(R.p).*10.^(2*(1:length(R.p)))];
p = [R.x0(1) R.x0(2) R.p];
if opts.MaxIter>0
    p = lsqnonlin(@uL2rdres,p,[],[],opts,R,u,L); % optimize
    R.x0 = p(1:2); R.x0 = R.x0(:);
    % R.p = sign(p(3:end)).*abs(p(3:end))./10.^(2*(1:length(p)-2));
    R.p = p(3:end);
end
if nargout>1 % corrected points
    for i=1:length(u)
        x{i} = ur2u(u{i},R); % undistorted points
    end
end
if nargout>3 % residuals
    for i=1:length(u)
        [e{i},z{i}] = uL2rdres(p,R,u(i),L(i));
    end
else
    for i=1:length(u)
        e{i} = uL2rdres(p,R,u(i),L(i));
    end    
end
