% [u,z] = X2u(X,C) - perspective projection
%
%  X      ... 3(4) x n space points, 3 -> 1's augmented
%  C      ... camera description
%  C.type ... camera type:
%
%   'P'   - C - 3 x 4 camera projection matrix  
%   'KRC' -
%                [1/bx s/f u0/f] - internal calibration matrix
%         C.K  = [0   1/by v0/f]   If f is in world units, true camera size
%                [0     0   1/f]   can be reconstruceted.
%         any non-zero multiple of C.K is equivalent for making projections
%
%         C.R  = 3x3 rotation matrix
%         C.C  = 3x1 camera center 
%   'KRCrp'
%         C.r  = 1xN [p1 p2 ...] polynomial distortion coefficients   
%   'KRCrd'
%         C.r  = 1x1 division model distortion parameter
%                model rd = (1-sqrt(1-4*C.r*ru^2))/(2*C.r*ru) 
%                      ru = rd * 1/(1+C.r*rd^2)
%                      ru - undisterted radius, rd - distorted radius
%   'OCV'
%         C.R = rotation          ... option 1  
%         C.C = 3x2 camera center ... option 1
%         C.A = [R -R*C;0 0 0 1]  ... option 2
%         C.r = radial distortion parameters (see OpenCV)
%         C.K = K matrix 
%
% Bakwards compatibility
%
%  C.type missing
%          C = P \in R^{3 x 4} or        
%  C.type missing
%          C.K  ... internal calibration matrix
%          C.E  ... camera euclidean coordinate system transformation
%          C.r  ... polynomial radial distortion
%
% u ... 3 x n image projections, affine coords. for finite points
% z ... z coordinate of the directiin vector (> for points infront of the
%                                             camera)

% (c) T.Pajdla, pajdal@cmp.felk.cvut.cz, 2015-07-12
function [u,z,d1,d2] = X2u(X,Ci)
% backwards compatibility of camera models
if ~isa(Ci,'struct')
    C.P = Ci; 
    C.type = 'P';
else
    C = Ci;
end
if ~isfield(C,'type')
    C.type = 'KRCrd';
    C.f = norm(C.E(3,1:3)); % focal length
    C.E = C.E/C.f; % normalize to get rotation 
    C.R = C.E(1:3,1:3); % rotation 
    C.C = -C.R'*C.E(:,4); % camera center
end
switch C.type
    case 'P'       
        cType = 'P';       rType = '';
    case 'KRC';
        cType = 'KRC';     rType = '';
    case 'KRCrp'   
        cType = 'KRC';     rType = 'rp';
    case 'KRCrd'   
        cType = 'KRC';     rType = 'rd';
    case 'CAHVOR'  
        cType = 'CAHVOR';  rType = '';
    case 'CAHVORe' 
        cType = 'CAHVORe'; rType = '';
    case 'OCV'
        cType = 'OCV';     rType = '';
end
% normalize format of X
if size(X,1)<4, 
    X = [X;ones(1,size(X,2))]; 
end
% projections
switch cType
    case 'P'
        u  = C.P*X;
        z = u(3,:);
        u  = u./(ones(3,1)*sqrt(sum(u.^2)));
        u3 = u(3,:);
        % normalize points at infinity to norm one
        u3(abs(u3)<1e-8) = 1;
        u  = u./(ones(3,1)*u3);
    case 'KRC'
        u = C.R*[eye(3) -C.C]*X;
        z = u(3,:); 
        u = u./(ones(3,1)*u(3,:));
        u = u(1:2,:);
        if any(abs(C.r)>10*eps)
            t = ones(1,size(u,2)); % the final radius parameter
            r = sqrt(sum(u.^2));
            for i=1:size(C.r,2)
                t = t + C.r(i)*r.^(2*i);
            end
            u = ([1;1]*t).*u;
        end
        u  = a2h(h2a(C.K*a2h(u)));
    case 'CAHVOR'
        if size(X,1)>3, X = X(1:3,:); end
        u = ones(3,size(X,2));
        oo = C.O*C.O'; % projector to unit O
        ox = eye(3)-oo; % perpediculator to unit O
        for i=1:size(X,2)
            Y = X(:,i)-C.C; % (2.51)
            t = sum((ox*Y).^2)/sum((oo*Y).^2); % (2.52)
            Z = Y+((t.^(1:size(C.R,2)))*C.R')*(ox*Y); % (2.53)
            z = Z(3,:);
            u(:,i) = [[C.H';C.V']*Z/(C.A'*Z);1]; % (2.54)
        end
    case 'CAHVORE'        
        if size(X,1)>3, X = X(1:3,:); end
        for i=1:size(X,2)
            % solve equation
            %
            % sin(t)*(b-(t/sin(t)-1)*(e1+e2*t^2+e3*t^4+...))-a*cos(t) = 0
            %
            % for t. By the Newton method. Initializaze: t = atan(a/b)
            %
            % a = ||(I-O*O')*(X-C)||
            % b = ||(O*O')*(X-C)|| = O'*(X-C)
            oo = C.O*C.O'; % projector to unit O
            ox = eye(3)-oo; % perpediculator to unit O 
            V = X(:,i)-C.C; 
            a = sqrt(sum((ox*V).^2)); 
            b = sqrt(sum((oo*V).^2));
            % solve for ty
            ty = atan2(a,b);
            for j=1:6
                f = ty.^(2*(0:size(C.E,1)-1))*C.E;
                g = ty/sin(ty)-1;
                e = sin(ty)*(b-f*g)-a*cos(ty);
                df = (2*(1:size(C.E,1)-1)).*(ty.^(2*(1:size(C.E,1)-1)-1));
                dg = (sin(ty)-ty*cos(ty))/(sin(ty)^2);
                de = cos(ty)*(b-f*g)-sin(ty)*(df*g+f*dg)+a*sin(ty);
                ty = ty-e/de;
            end
            f = ty.^(2*(0:size(C.E,1)-1))*C.E;
            s  = (ty/sin(ty)-1)*f; % (2.124)
            Cp = C.C + s*C.O; % (2.123)
            Y = X(:,i)-Cp; % (2.125)
            xi = sin(C.L*t)/(L*cos(max(0,L*t))); % (2.127)
            ttz = (1+(t.^(2*(0:size(C.R,2)))*C.R))*xi; % (2.126)
            tty = tan(ty);
            Z = ((ttz/tty)*ox+oo)*Y; % (2.129)
            z = Z(3,:);
            u(:,i) = [C.H';C.V']*Z/(C.A'*Z); % (2.54)
        end
    case 'OCV'
        if ~isfield(C,'A')
            C.A = [C.R -C.R*C.C;[0 0 0 1]]; % Metric projection matrix
        end
        Z = C.A*X;
        z = Z(3,:);
        p = h2a(h2a(Z)); % projection to the metric image
        r2 = sum(p.^2); % squared radius
        d1 = 1 + C.r(1)*r2 + C.r(2)*r2.*r2 + C.r(5)*r2.*r2.*r2;
        if (numel(C.r) == 8)
            d2 = 1 + C.r(6)*r2 + C.r(7)*r2.*r2 + C.r(8)*r2.*r2.*r2;
        else
            d2 = 1;
        end
        % radial distortion multiplier
        pp(1,:) = (d1./d2).*p(1,:)+2*C.r(3)*p(1,:).* p(2,:)+ C.r(4).*(r2 + 2*p(1,:).*p(1,:));
        pp(2,:) = (d1./d2).*p(2,:)+2*C.r(4)*p(1,:).* p(2,:)+ C.r(3).*(r2 + 2*p(2,:).*p(2,:));
        % to the image coordinate system
        u(1,:) = C.K(1,1)*pp(1,:) + C.K(1,3);
        u(2,:) = C.K(2,2)*pp(2,:) + C.K(2,3);
end