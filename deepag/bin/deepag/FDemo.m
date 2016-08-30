% deepag.m - Deep Learning 4 Algebraic Geometry - FDemo, nogui
% T. Pajdla, pajdla@gmail.cz
% 1 Jan 2016 - 31 Dec 2016
% 
% 
function [fnb,fgnb] = FDemo()
    %% Initialize
    deepagpaths;          % set paths to libraries
    ps.Data    = 'paris'; % Data path
    ps.FDemo   = true; % Cameras & F computation demo
    ps.ImNS    = 1.0;   % Image noise ~ N(0,ps.ImNS)
    ps.plot    = false;   % Image noise ~ N(0,ps.ImNS)
    deepagInit
    %%
    if ps.FDemo
        % focal lengths
        f1 = 900;
        f2 = 1100;
        sc = max(f1,f2);
        % 3D points
        X = 2*sc*[0 1 1 0 0 1 2 0
                 0 0 1 1 0 0 1 1
                 0 0 0 0 1 2 1 2]+repmat([-sc;-sc;2*sc],1,8);
        % Internal camera calibration     
        K1 = [f1  0    0 
              0   f1   0 
              0   0    1];
        K2 = [f2 0  0
              0  f2 0
              0  0  1];
        % Rotations  
        R1 = eye(3);
        R2 = a2r(rand(3,1),pi/10);
        % Projection matrices
        P1 = K1*R1*[eye(3)     [-sc;0;0]];
        P2 = K2*R2*[eye(3)  sc*[rand(2,1);0]];
        % Image points
        u1 = X2u(X,P1);
        u2 = X2u(X,P2);
        % add image noise
        u1(1:2,:)=u1(1:2,:)+ps.ImNS*randn(size(u1(1:2,:)));
        u2(1:2,:)=u2(1:2,:)+ps.ImNS*randn(size(u2(1:2,:)));
        % Fundamental matrix
        [F,A]  = uu2F({u1,u2},{'[-1,1]','Free'}); % from image measurements
        F = F/norm(F);
        Fg = PP2F(P1,P2); % from projection matrices, ground truth
        Fg = Fg/norm(Fg);
        if ps.plot
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
        end
        % Compute the difference in normalized image coordinates
        Fn = inv(A{2})'*F*inv(A{1});
        Fn = Fn/norm(Fn);
        Fgn = inv(A{2})'*Fg*inv(A{1});
        Fgn = Fgn/norm(Fgn);
        if ps.plot
            % Plot the difference
            subfig(3,4,6);
            e = [Fgn(:)/norm(Fgn)-Fn(:)/norm(Fn)];
            plot(e,'.-r');axis tight;xlabel('Fn(:)');ylabel('error');title('Normalized: F(:)-PP2F(P1,P2)(:)');
        end
        % Form the feature vector and the ground truth feature vector
        f = Fn(:)/norm(F,'fro'); [~,mi] = max(abs(f)); f = f*sign(f(mi));
        fg = Fgn(:)/norm(Fg,'fro'); [~,mi] = max(abs(fg)); fg = fg*sign(fg(mi));
        if ps.plot
            subfig(3,4,7);
            plot(1:length(f),f+1,'.-',1:length(fg),fg+1,'.-g');axis tight;xlabel('element');ylabel('value');title('b - f, g - fg');
        end
        % Bougnoux formula aplied on normalized data and recomputed to original data
        ff = F2f1f2(F);
        fn = F2f1f2(Fn);
        fgn = F2f1f2(Fgn);
        fnb = fn*diag([1/A{1}(1) 1/A{2}(1)]);
        fgnb = fgn*diag([1/A{1}(1) 1/A{2}(1)]);
        if ps.plot
            % plot estimated f's
            plot([1 2],[[f1 f2];ff;fnb]','.','markersize',20);
            xlabel('focal length index'); ylabel('pixels');
            title(sprintf('b-truth, g-est, r-norm im(f) = %.1f %.1f',imag(fnb(1)),imag(fnb(2))));
        end
    end
end