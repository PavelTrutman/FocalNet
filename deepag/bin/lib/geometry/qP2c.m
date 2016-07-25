% [c,A] = qP2c(u,K) - the projection of the center of a circle observed as an ellise in an mage with known K matrix 
%
% u = 2 x n image coordinates of point on the observed ellipse
% K = camera internal K-matrix
% c = 2 x 1 projection of the center of the observed circle
% A = rectification matrix from the ellipse to the circle 
function [cc,A,f] = qP2c(u,K,dbg)

if nargin>0
    if nargin<3, dbg = false; end
   [p,q] = fitellipse(u); % fit ellipse
    c = p(1:2)'; % get its center
    v = K \ a2h(c); % the vector approximating the direction to the circle center
    v = v/vnorm(v); % unit
    R1 = xy2R3(v,[0;0;1]); % camera rotation bringing the ellipse into the image center
    u1 = h2a(R1*inv(K)*a2h(u)); % points in the center of the image
    [p1,q1] = fitellipse(u1); % fit ellipse to get the ellipse matrix
    [p1,A1] = q2par(q1); % get the rectification transform A2
    u2 = h2a(inv(A1)*R1*inv(K)*a2h(u)); % points rectified to have the shorter semiaxis aligned with the y direction in the image    
    [p2,q2] = fitellipse(u2); % fit ellipse to get the ellipse matrix  
    a2 = acos(min(p2(3:4))/max(p2(3:4)));
    A21 = a2r([1;0;0],a2);
    A22 = a2r([1;0;0],-a2);
    u31 = h2a(A21*inv(A1)*R1*inv(K)*a2h(u)); % points projected onto the image plene roughly parallel with the plane of the target 
    u32 = h2a(A22*inv(A1)*R1*inv(K)*a2h(u)); % points projected onto the image plene roughly parallel with the plane of the target 
    [p31,q31] = fitellipse(u31); % fit ellipse to get the ellipse matrix     
    [p32,q32] = fitellipse(u32); % fit ellipse to get the ellipse matrix
    x31 = h2a(inv(A21*inv(A1)*R1*inv(K))*a2h(p31(1:2)'));
    x32 = h2a(inv(A22*inv(A1)*R1*inv(K))*a2h(p32(1:2)'));
    [xd,xdi] = min(vnorm([c c]-[x31 x32]));
    if xdi==1
        A = A21*inv(A1)*R1*inv(K);
        cc = x31;
    else
        A = A22*inv(A1)*R1*inv(K);
        cc = x32;
    end
    if dbg 
        f(1)=subfig(2,3,1); 
        plot3d(u,'.'); axis equal; hold
        plot3d(h2a(sampleq(q,2*size(u,2))),'r');
        plot3d(c,'.r');       
        f(2)=subfig(2,3,2);
        plot3d(u1,'b.');hold; axis equal
        plot3d(h2a(sampleq(q1,2*size(u,2))),'b');        
        plot3d(h2a(R1*inv(K)*a2h(c)),'.b');
        plot3d(p1(1:2),'+b');
        plot3d(u2,'.m');
        plot3d(h2a(sampleq(q2,2*size(u,2))),'m');        
        plot3d(h2a(inv(A1)*R1*inv(K)*a2h(c)),'.m');
        plot3d(p2(1:2)','+m');
        plot3d(u31,'.k');
        plot3d(h2a(sampleq(q31,2*size(u,2))),'k');        
        plot3d(h2a(A21*inv(A1)*R1*inv(K)*a2h(c)),'.k');
        plot3d(p31(1:2)','+k');
        plot3d(u32,'.k');
        plot3d(h2a(sampleq(q32,2*size(u,2))),'k');        
        plot3d(h2a(A22*inv(A1)*R1*inv(K)*a2h(c)),'.k');
        plot3d(p32(1:2)','xk');        
        figure(f(1))
        plot3d(x31,'+k');
        plot3d(x32,'xk');
    end
else % demo
    N = 30; % the number of points on the circle
    X = [h2a(sampleq(diag([1 1 -10^2]),N));zeros(1,N)];
    K = [1000 0 500;0 1000 500;0 0 1]; % K matrix
    R = a2r([1;0;0],-pi/4.1); % rotation
    C = [0;-300;300];
    P = K*R*[eye(3) -C]; % projection matrix
    u = h2a(P*a2h(X)); % project
    [c,A,f] = qP2c(u,K,true);
    figure(f(1));plot3d(h2a(P*a2h([0;0;0])),'og');plot3d(c,'.g'); plot([0 1000 1000 0 0],[0 0 1000 1000 0]); 
    figure(f(2));plot3d(h2a(A*P*a2h([0;0;0])),'og');
    vnorm(c-h2a(P*a2h([0;0;0])))
end






