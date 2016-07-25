% e = ux2rdres(p,R,u,X) - radial distortion estimation from homographies - residuals
% 
% R     ... radial distortion model, see ur2u 
% p     ... [c(1) c(2) rp] ... model parameters
%            c = center
%            rp = distortion parameters, see ur2u
% u     ... u{i} = 2 x n projections of the i-th target 
% X     ... X{i} = 2 x n coordinates of the i-th target 
% e     ... homography fit error
% z     ... 4 x n [the closest points z on the line; x]
%
% Assumes that x and u were normalized.
%
% See also ux2rdc

% (c) T.Pajdla, cmp.felk.cvut.cz, May 1, 2006
function [e,z] = ux2rdres(p,R,u,X)

%% the standard radial distortion model, see u2rd
R.x0(1:2) = p(1:2);
R.p = p(3:end);
z = cell(1,length(u));
%% homography fit residuals
for i=1:length(u) % target projections
    x = a2h(X{i});
    y = a2h(ur2u(u{i},R)); % undistorted points
    [H,e{i}]=xy2H(x,y,{'DLT','y-y'}); % best homography fit
    z{i} = h2a(H*x);
end
%% reshape the residuals into a vector
e = cat(2,e{:});
e = e(:);
if nargout>1
    z = cat(2,z{:});
end
