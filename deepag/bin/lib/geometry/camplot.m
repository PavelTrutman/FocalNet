% h = camplot(P[,f,c,u,t]) - perspective camera plotting
%
% P = Camera structure P.K, P.R, P.C see X2u.m
%     or 3x4 perspective camea projection matrix P = K*R*[I|-C] with rotation R, pose C and
%
%          [1/bx s/f  u0/f] 
%     K  = [0    1/by v0/f]   
%          [0     0    1/f]   
%
%     determining the camera coordinate system \beta. 
%
% f = overrides focal length in K by f, if present and not empty
% c = image corners in \alpha, uses [0 2*u0 2*u0 0;0 0 2*v0 2*v0] if missing
% 
% h = plot handles 
%     h(1)      = image frame
%     h(2:5)    = frame to center connectors 
%     h(6)      = center to principal point connector
%     h(7)      = center dot
%     h(8)      = center circle
%     h(9)      = first axis
%     h(10)     = second axis
%     h(11:end) = image points

% Tomas Pajdla (pajdla@cmp.felk.cvut.cz)
% 2015-09-04
function h = camplot(P,f,c,u,t)
if isstruct(P)
    K = P.K;
    R = P.R;
    C = P.C;
    P = KRC2P(P);
else
    [K,R,C] = P2KRC(P); % camera parameters
end
if all(isfinite(P(:)))
    if nargin>1 % override f
        if ~isempty(f)
            K = K/K(3,3)/f;
        end
    end
    o_a = [0;0]; % origin of \aplha in \alpha
    p_a  = K(1:2,3)/K(3,3); % principal point in \alpha
    p_d  = [R' C]*a2h(K\a2h(p_a)); % principal point in \delta
    if nargin>2 && ~isempty(c)
        x_a = c(:,[1 2 3 4 1]);
        s = vnorm([x_a(:,2)-x_a(:,1) x_a(:,4)-x_a(:,1)])/2; 
        % scaled basic vectors of \alpha in \delta 
        a_d = [R' C]*a2h(K\a2h([min(s) 0 
                                0      min(s)]));
    else
        x_a = [0 2*p_a(1) 2*p_a(1) 0        0
               0 0        2*p_a(2) 2*p_a(2) 0]; % outer frame in \alpha
        % scaled basic vectors of \alpha in \delta    
        a_d = [R' C]*a2h(K\a2h([min(p_a) 0 
                                0      min(p_a)]));
    end
    o_d = [R' C]*a2h(K\a2h(o_a)); % origin of \aplha in \delta
    x_d  = [R' C]*a2h(K\a2h(x_a)); % outer frame in \delta

    h = plot3d(x_d,'k');
    if ~ishold, hold; unhold=true; else unhold=false; end
    h = [h plot3d([C x_d(:,1)],'k')];
    h = [h plot3d([C x_d(:,2)],'k')];
    h = [h plot3d([C x_d(:,3)],'k')];
    h = [h plot3d([C x_d(:,4)],'k')];    
    h = [h plot3d([C p_d],'-g')];
    h = [h plot3d(C,'g.','markersize',15)];
    h = [h plot3d(C,'ok')];
    h = [h plot3d([o_d a_d(:,1)],'b','linewidth',2)];
    h = [h plot3d([o_d a_d(:,2)],'r','linewidth',2)];
    if nargin>3
        u = [R' C]*a2h(K\a2h(u)); % ray direction vectors
        h = [h plot3d(u,'.b')];
    end
    if unhold, hold; end
else
    h = []; % to work in cellfun for non-valid cameras filled with NaNs
end