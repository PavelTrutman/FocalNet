% [e, z] = uL2rdres(p,R,u,L) - radial distortion estimation from lines - residuals
% 
% R     ... radial distortion model, see ur2u 
% p     ... [c(1) c(2) p] ... model parameters, see ur2u
% u     ... u{i} = 2 x n projections of the i-th target 
% L     ... L{i}{j} = indices of points lying on a stright line in the tartget
%       ... normalization of ccordinates, see R.A in ur2u 
% e     ... line fit errors
% z     ... 4 x n [the closest points z on the line; x]
% 
% Assumes that x and u were normalized.
%
% See also ux2rdc

% (c) T.Pajdla, cmp.felk.cvut.cz, May 1, 2006
function [e,z] = uL2rdres(p,R,u,L)

%% the standard radial distortion model, see u2rd
R.x0(1:2) = p(1:2);
R.p = p(3:end);
%% line fit residuals
k = 1; 
for i=1:length(u) % target projections
    x = ur2u(u{i},R); % undistorted points
    for j=1:length(L{i}) % lines
        y = x(:,L{i}{j}); % points on the line j
        if nargout>1
            [l,e{k},z{k}] = l2Dfit(y); % fit a line
        else
            [l,e{k}] = l2Dfit(y); % fit a line
        end
        if 0
            v = l(1:2)*e{k};
            figure; plot(y(1,:),y(2,:),'.'); hold; plotline(l); 
            quiver(y(1,:),y(2,:),v(1,:),v(2,:),1);
            title(num2str(mean(abs(e{k}))));
            axis equal
            pause
        end
        k = k + 1;
    end
end
%% reshape the residuals into a vector
e = cat(2,e{:});
e = e(:);
if nargout>1
    z = cat(2,z{:});
end
