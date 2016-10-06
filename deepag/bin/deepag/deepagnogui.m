% deepag.m - Deep Learning 4 Algebraic Geometry - nogui
% T. Pajdla, pajdla@gmail.cz
% 1 Jan 2016 - 31 Dec 2016
%
%% Initialize
deepagpaths;          % set paths to libraries
ps.Data    = 'paris'; % Data path
ps.FDemo   = false;    % Cameras & F computation demo
ps.SceneGen= true;   % Scene generator
ps.ptN     = 7;       % The number of points for F computation
ps.ImNS    = 0.0;     % Image noise ~ N(0,ps.ImNS)
deepagInit
%%
if ps.FDemo
   % focal lengths
    f1 = 900;
    f2 = 1100;
    sc = max(f1,f2);
    % 3D points
    X = [0 1 1 0 0 1 2 1 
         0 0 1 1 0 0 1 0
         0 0 0 0 1 2 1 1];
    X = 2*sc*(X+rand(size(X)))+repmat([-sc;-sc;2*sc],1,size(X,2));
    % Internal camera calibration     
    K1 = [f1  0    0 
          0   f1   0 
          0   0    1];
    K2 = [f2 0  0
          0  f2 0
          0  0  1];
    % Rotations  
    R1 = a2r(rand(3,1),pi/8);
    R2 = a2r(rand(3,1),pi/8);
    % Projection matrices
    P1 = K1*R1*[eye(3)     [-10*sc;0;0]];
    P2 = K2*R2*[eye(3)  10*sc*([rand(2,1);0])];
    % Image points
    u1 = X2u(X,P1);
    u2 = X2u(X,P2);
    % add image noise
    u1(1:2,:)=u1(1:2,:)+ps.ImNS*randn(size(u1(1:2,:)));
    u2(1:2,:)=u2(1:2,:)+ps.ImNS*randn(size(u2(1:2,:)));
    % Fundamental matrix
    clear e
    [F,A]  = uu2F({u1,u2},{'None','|F|=0'}); % from image measurements
    for i=1:size(F,3)
        F(:,:,i) = F(:,:,i)/norm(F(:,:,i));
        e(i) = max(abs(sum(u2.*(F(:,:,i)*u1))));
    end    
    [~,ie] = min(e);
    F = F(:,:,ie);    
    Fg = PP2F(P1,P2); % from projection matrices
    Fg = Fg/norm(Fg);
    % Plot cameras & 3D points
    subfig(3,4,1);
    camplot(P1,[],f1/2*[-1 1 1 -1;-1 -1 1 1]);hold; camplot(P2,[],f2/2*[-1 1 1 -1;-1 -1 1 1]); 
    plot3d(X,'.');
    axis equal
    title('Cameras & 3D points (''.b'')')
    % Plot images
    subfig(3,4,2);plot3d(u1(1:2,:),'.');axis image;title('Image 1'); 
    subfig(3,4,3);plot3d(u2(1:2,:),'.');axis image;title('Image 2');
    % Plot algebraic error
    subfig(3,4,4);
    e = sum(u2.*(F*u1));
    plot(e,'.-r');axis tight;xlabel('point #');ylabel('error');title('Algebraic error e = u2''*F*u1''');
    % Compare the computed F with the ground truth Fg
    % In image coordinates
    subfig(3,4,5);
    e = [Fg(:)/Fg(3,3)-F(:)/F(3,3)];
    plot(e,'.-r');axis tight;xlabel('F(:)');ylabel('error');title('F(:)-PP2F(P1,P2)(:)');    
    % Compute the dinfference in normalized image coordinates
    Fn = inv(A{2})'*F*inv(A{1});
    Fn = Fn/norm(Fn);
    Fgn = inv(A{2})'*Fg*inv(A{1});
    Fgn = Fgn/norm(Fgn);
    % Plot the difference
    subfig(3,4,6);
    e = [Fgn(:)/norm(Fgn)-Fn(:)/norm(Fn)];
    plot(e,'.-r');axis tight;xlabel('Fn(:)');ylabel('error');title('Normalized: F(:)-PP2F(P1,P2)(:)');
    % Form the feature vector and the ground truth feature vector
    f = Fn(:)/norm(F,'fro'); [~,mi] = max(abs(f)); f = f*sign(f(mi));
    fg = Fgn(:)/norm(Fg,'fro'); [~,mi] = max(abs(fg)); fg = fg*sign(fg(mi));
    subfig(3,4,7);
    plot(1:length(f),f+1,'.-',1:length(fg),fg+1,'.-g');axis tight;xlabel('element');ylabel('value');title('b - f, g - fg');
    % Bougnoux formula aplied on normalized data and recomputed to original data
    ff = F2f1f2(F)
    fn = F2f1f2(Fn);
    fgn = F2f1f2(Fgn);
    fnb = fn*diag([1/A{1}(1) 1/A{2}(1)])
    fgnb = fgn*diag([1/A{1}(1) 1/A{2}(1)]);
    % plot estimated f's
    plot([1 2],[[f1 f2];ff;fnb]','.','markersize',20);
    xlabel('focal length index'); ylabel('pixels');
    title(sprintf('b-truth, g-est, r-norm im(f) = %.1f %.1f',imag(fnb(1)),imag(fnb(2))));
end
if ps.SceneGen
    sceneType = {'randombox' 'random'};
    pixel = 1/1000;
    noise = 1*pixel;
    Npoints = 20;
    Ncams = 2;
    samplesize = 8;
    f1 = 900*pixel;
    f2 = 1100*pixel;
    Kgt1 = diag([f1 f1 1]);
    Kgt2 = diag([f2 f2 1]);
    gtk1 = 0;
    gtk2 = 0;
    [Pgt M m mgt] = GenerateScene(Npoints, [4000*pixel 4000*pixel], Ncams, 5000*pixel, 8000*pixel, 0, noise, [Kgt1;Kgt2], sceneType, [], [], [gtk1,gtk2], true);
    Kgt{1} = Kgt1; Kgt{2}=Kgt2;
    %ShowCameras(Pgt, Kgt, m, M, true, false, true, 1:7, mgt);
    sample=randperm(size(m{1},2),samplesize);
    u={m{1}(:,sample)/pixel m{2}(:,sample)/pixel};
    [F,A]=F_features(u{1}, u{2});
    estion=F2f1f2(reshape(F,3,3));
    estion=estion*diag([1/A{1}(1) 1/A{2}(1)])
end
