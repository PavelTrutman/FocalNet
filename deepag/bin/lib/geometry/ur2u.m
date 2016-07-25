% [x,r,c,K] = ur2u(u,R) - radial distorion correction 
% 
% u ... 2 x n image points
% R ... radial distortion descriptor
%
%       if a structure:
%
%       R.x0  ... 2 x 1 distortion center ([0;0] implicit if field r.x0 missing)
%       R.p   ... 1 x k ... k radial distortion coefficients
%       R.m   ... model type
%                 1 ... polynomial (implicit if field r.m missing)
%                       c = 1 + p(1)*r.^2 + p(2)*r.^4 + ...
%                       x = inv(R.A)*(((ones(size(x0,1),1)*c) .* R.A*(u-(R.x0*ones(1,size(u,2)))))) + (R.x0*ones(1,size(u,2)));
%                       c ... 1 x size(u,2) coefficients
%                       r ... 1 x size(u,2) radii
%                 2 ... the division model
%                       r = r' * 1/(1+a*r'^2) reverse by r' = (1-sqrt(1-4*a*r^2))/(2*a*r)
%                       c ... 1 x size(u,2) coefficients
%                       r ... 1 x size(u,2) radii
%                 3 ... OpenCV model 
%                       R.fm  ... undistorted image focal length multiplicationi constant (defaul = 1)
%                       R.ims ... undistorted image image scale
%                       multiplication consant (default = 1)
%                       r = c = [];
%                       K ... output K matrix
%       R.A   ... 2 x 2 matrix, coordinate normalization of a vector 
%
% x ... corrected points

% (c) T.Pajdla, www.neovision.cz, Aug 4, 2005
function [u,r,c,Ko] = ur2u(u,R)

% the implicit model 
if ~isfield(R,'m')
    if isfield(R,'type') % get parameters form the camera descripton used in X2u
        switch R.type
            case 'KRCrd', 
                R.m = 2; % divison model
                R.x0 = R.K(1:2,3); % distortion center is in the main point
            case 'KRCrp', 
                R.m = 1;
                R.x0 = R.K(1:2,3); % distortion center is in the main point
            case 'OCV'
                R.m = 3;
            otherwise
                error('not implemented');
        end         
    else
        R.m = 1;
    end
end
% implicit center
if ~isfield(R,'x0')
    R.x0 = [0;0];
else
    R.x0 = R.x0(:);
end
% implicit normalization
if ~isfield(R,'A')
    R.A = eye(2);
end
% undistort
o1 = ones(1,size(u,2));
switch R.m
    case 1 % polynomial model        
        x0 = R.x0*o1;
        u  = R.A*(u-x0);
        r  = sqrt(sum(u.^2));
        c  = 0*o1; % zeros
        for i=1:length(R.p)
            c = c + R.p(i)*r.^(2*i);
        end
        c = c + 1;
        u = inv(R.A)*((ones(size(x0,1),1)*c).*u)+x0;
    case 2 % the division model
        x0 = R.x0*o1;
        u  = R.A*(u-x0);
        r  = sqrt(sum(u.^2));
        c  = 0*o1; % zeros
        for i=1:length(R.p)
            c = c + R.p(i)*r.^(i+1);
        end        
        c = 1./(c + 1);
        u = inv(R.A)*((ones(size(x0,1),1)*c).*u)+x0;
    case 3 % OCV model
        if ~isfield(R,'ims'), R.ims = 1; end
        if ~isfield(R,'fm'), R.fm = 1; end
        dm = (1/max(abs(R.K([1 5])))/1000)^2; % a fraction of pixel/f squared 
        K = R.K; % input K
        Ko = K; % output K
        Ko(1,1) = Ko(1,1)*R.fm; 
        Ko(2,2) = Ko(2,2)*R.fm;
        Ko(1:2,1:3) = Ko(1:2,1:3)*R.ims;
        dist = R.r;
        u_dist = u;
        u = zeros(size(u_dist));
        x = u; d = u;
        ps = [0;0];
        for i = 1:size(u,2)
            ps(1) = (u_dist(1,i) - K(1,3)) / K(1,1); % into the camera coordinate system
            ps(2) = (u_dist(2,i) - K(2,3)) / K(2,2);
            p = ps;
            for j = 1:100 % itteratively solve for the inverse of distortion
                r2 = p(1) * p(1) + p(2) * p(2);
                if (numel(dist) == 8)
                    icdist = (1.0 + ((dist(8)*r2 + dist(7))*r2 + dist(6))*r2)/(1.0 + ((dist(5)*r2 + dist(2))*r2 + dist(1))*r2);
                else
                    icdist = 1.0/(1 + ((dist(5)*r2 + dist(2))*r2 + dist(1))*r2);
                end
                dx = 2 * dist(3) * p(1) * p(2) + dist(4) * (r2 + 2 * p(1) * p(1));
                dy = 2 * dist(4) * p(1) * p(2) + dist(3) * (r2 + 2 * p(2) * p(2));
                po = p;
                p(1) = (ps(1) - dx) * icdist;
                p(2) = (ps(2) - dy) * icdist;
                if (p(1)-po(1))^2+(p(2)-po(2))^2<dm, break; end % exitif the change is small
            end
            n(i) = j;
            d(:,i) = p-po;
            x(:,i) = p;
            p(1) = p(1) * Ko(1,1) + Ko(1,3);
            p(2) = p(2) * Ko(2,2) + Ko(2,3);
            u(:,i) = p;
        end
        c = []; r = [];
end

