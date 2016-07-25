% [X,v,z,e,a] = uu2X(u,P,m) - optimal point reconstruction from projections
%
% u = image matches [[u1;v1;...;u2;v2;...] ...] 
% P = projection matrices [P1;P2;...]
% X = optimal reconstruction of points minimizing the sum of squared reprojection errors
% v = optimal projection of X
% z = 2 x n z coordinates of projections in the camera coordinate system (z>0 is infront of cameras)
% e = reprojection error ||v-u||
% a = apical angles for multi-view reconstruction
% m = method 
%     two-view
%     'PROJ' = optimasl projective reconstruction [HZ-2003] (used if missing)
%     'TRAN' = shortest transversal
%     multi-view
%     two view reconstruction for the rays with the apical angle at X colsest to pi/3

% pajdla@cmp.felk.cvut.cz, 2015-09-05
function [X,v,z,e,a] = uu2X(u,P,method)
if nargin>0
    if size(u,1)>4 % multiview largest angle
        X = zeros(3,size(u,2)); % 3D points
        a = zeros(1,size(u,2)); % maximal apical angles
        ix = nchoosek(1:size(P,1)/3,2); % all pairs
        for i=1:size(u,2)
            Y = zeros(3,size(ix,1)); % points for all pairs
            ai = zeros(1,size(ix,1)); % apical angles for all pairs
            for j=1:size(ix,1)
                Pj = [P(3*(ix(j,1)-1)+(1:3),:); P(3*(ix(j,2)-1)+(1:3),:)];
                uij = [u(2*(ix(j,1)-1)+(1:2),i); u(2*(ix(j,2)-1)+(1:2),i)];
                Y(:,j) = uu2X(uij,Pj,'TRAN'); % reconstruct
                aa = PPX2ae(Pj(1:3,:),Pj(4:6,:),Y(:,j)); % apical angle
                ai(j) = aa(3); % at X
            end
            [~,mi] = min(abs(pi/3-ai));
            a(i) = ai(mi);
            X(:,i) = Y(:,mi);
        end
        v = zeros(size(u));
        z = zeros(size(P,1)/3,size(u,2));
        e = zeros(size(P,1)/3,size(u,2));
        for i=1:size(P,1)/3
            [uu,z(i,:)] = X2u(X,P(3*(i-1)+(1:3),:));
            v(2*(i-1)+(1:2),:) = h2a(uu);
            e(i,:) = vnorm(v(2*(i-1)+(1:2),:)-u(2*(i-1)+(1:2),:));
        end
    else % two views
        if nargin<3
            method = 'PROJ';
        end
        P1 = P(1:3,:);
        P2 = P(4:6,:);
        % remove scale & offset in 3D coordinate system
        P1 = P1/vnorm(P1(3,1:3));
        P2 = P2/vnorm(P2(3,1:3));
        sP = mean([P1(1:3,1:3)\P1(:,4) P2(1:3,1:3)\P2(:,4)],2);
        TP = [eye(3) -sP; 0 0 0 1];
        P1 = P1*TP;
        P2 = P2*TP;
        switch method
            case 'PROJ' % optimal projective
                % prepare
                A0 = [0 0 0 0 0 0;1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0];
                B0 = eye(6);
                %
                F = PP2F(P1,P2);
                F0 = F;
                v = zeros(size(u)); % corrected points
                X = zeros(3,size(u,2)); % 3D points
                for i=1:size(u,2) % for all correspondences
                    % algorithm [HZ-2003 p.318]
                    u1 = u(1:2,i);
                    u2 = u(3:4,i);
                    % change the coordinate system
                    iT1  = [1 0 u1(1)
                        0 1 u1(2)
                        0 0    1];
                    iT2 =  [1 0 u2(1)
                        0 1 u2(2)
                        0 0    1];
                    F  = iT2'*F0*iT1;
                    F  = F/norm(F);
                    e1 = null3x3r2(F); %e1 = e1/vnorm(e1(1:2));
                    e2 = null3x3r2(F');%e2 = e2/vnorm(e2(1:2));
                    R1 = [e1(1) e1(2) 0;
                        -e1(2) e1(1) 0;
                        0     0  1];
                    R2 = [e2(1) e2(2) 0;
                        -e2(2) e2(1) 0;
                        0     0  1];
                    F  = R2*F*R1';
                    % get parameters
                    a = F(2,2);
                    b = F(2,3);
                    c = F(3,2);
                    d = F(3,3);
                    f1 = e1(3);
                    f2 = e2(3);
                    % coefficients of the polynomial - slow
                    % p6 = -a^2*c*d*f1^4+a*b*c^2*f1^4;
                    % p5 = -a^2*d^2*f1^4+b^2*c^2*f1^4+c^4*f2^4+2*a^2*c^2*f2^2+a^4;
                    % p4 = -a*b*d^2*f1^4+b^2*c*d*f1^4+4*c^3*d*f2^4-2*a^2*c*d*f1^2+4*a^2*c*d*f2^2+2*a*b*c^2*f1^2+4*a*b*c^2*f2^2+4*a^3*b;
                    % p3 = 6*c^2*d^2*f2^4-2*a^2*d^2*f1^2+2*a^2*d^2*f2^2+8*a*b*c*d*f2^2+2*b^2*c^2*f1^2+2*b^2*c^2*f2^2+6*a^2*b^2;
                    % p2 = 4*c*d^3*f2^4-2*a*b*d^2*f1^2+4*a*b*d^2*f2^2+2*b^2*c*d*f1^2+4*b^2*c*d*f2^2-a^2*c*d+4*a*b^3+a*b*c^2;
                    % p1 = d^4*f2^4+2*b^2*d^2*f2^2-a^2*d^2+b^4+b^2*c^2;
                    % p0 = -a*b*d^2+b^2*c*d;
                    a2=a^2;b2=b^2;c2=c^2;d2=d^2;f12=f1^2;f14=f12^2;f22=f2^2;f24=f22^2;
                    p6 = -a2*c*d*f14+a*b*c2*f14;
                    p5 = -a2*d2*f14+b2*c2*f14+c2^2*f24+2*a2*c2*f22+a2^2;
                    p4 = -a*b*d2*f14+b2*c*d*f14+4*c^3*d*f24-2*a2*c*d*f12+4*a2*c*d*f22+2*a*b*c2*f12+4*a*b*c2*f22+4*a^3*b;
                    p3 = 6*c2*d2*f24-2*a2*d2*f12+2*a2*d2*f22+8*a*b*c*d*f22+2*b2*c2*f12+2*b2*c2*f22+6*a2*b2;
                    p2 = 4*c*d^3*f24-2*a*b*d2*f12+4*a*b*d2*f22+2*b2*c*d*f12+4*b2*c*d*f22-a2*c*d+4*a*b^3+a*b*c2;
                    p1 = d2^2*f24+2*b2*d2*f22-a2*d2+b2^2+b2*c2;
                    p0 = -a*b*d2+b2*c*d;
                    % solve by companion matrix
                    %A = [0 0 0 0 0 -p0
                    %     1 0 0 0 0 -p1
                    %     0 1 0 0 0 -p2
                    %     0 0 1 0 0 -p3
                    %     0 0 0 1 0 -p4
                    %     0 0 0 0 1 -p5];
                    A = A0; A(:,6) = -[p0;p1;p2;p3;p4;p5];
                    %B = [1 0 0 0 0 0
                    %     0 1 0 0 0 0
                    %     0 0 1 0 0 0
                    %     0 0 0 1 0 0
                    %     0 0 0 0 1 0
                    %     0 0 0 0 0 p6];
                    B = B0; B(6,6) = p6;
                    t = eig(A,B);
                    % select real roots & add the t->infty
                    ix = abs(imag(t))<eps('single');
                    t = t(ix);
                    if abs(f1)>eps('single')
                        t = [t; 1/f1^2+c^2/(a^2+f2^2*c^2)];
                    end
                    % find the minimizer
                    e = t.^2./(1+f1^2*t.^2) + (c*t+d).^2./((a*t+b).^2+f2^2*(c*t+d).^2);
                    [~,im] = min(e);
                    t = t(im);
                    % get the optimal projections
                    l1 = [t*f1;1;-t]; % epipolar line in image 1
                    l2 = [-f2*(c*t+d)';a*t+b;c*t+d]; % epipolar line in image 2
                    v1 = [-l1(1)*l1(3);-l1(2)*l1(3);l1(1)^2+l1(2)^2]; % closest point to the origin
                    v2 = [-l2(1)*l2(3);-l2(2)*l2(3);l2(1)^2+l2(2)^2]; % closest point to the origin
                    % to the original coordinate system
                    v1 = h2a(iT1*R1'*v1);
                    v2 = h2a(iT2*R2'*v2);
                    v(:,i) = [v1;v2];
                    % 3D points
                    if true
                        % A = [xx(a2h(v1))*P1; xx(a2h(v2))*P2]; % slow
                        A = [[0 -1 v1(2);1 0 -v1(1);-v1(2) v1(1) 0]*P1
                            [0 -1 v2(2);1 0 -v2(1);-v2(2) v2(1) 0]*P2];
                        if all(isfinite(A(:)))
                            [~,~,x] = svd(A,'econ');
                            X(:,i) = x(1:3,end)/x(4,end);
                        else % no real finite solution
                            x = nan(4,1);
                            X(:,i) = nan(3,1);
                        end
                    else
                        C1 = -inv(P1(1:3,1:3))*P1(:,4);
                        C2 = -inv(P2(1:3,1:3))*P2(:,4);
                        A = [inv(P1(1:3,1:3))*a2h(v1) -inv(P2(1:3,1:3))*a2h(v2)];
                        b = C2-C1;
                        a = A\b;
                        x = a(1)*a2h(v1);
                        X(:,i) = x;
                    end
                end
            case 'TRAN'
                C1 = -P1(:,1:3)\P1(:,4); % centers
                C2 = -P2(:,1:3)\P2(:,4);
                T = C2-C1;
                x1 = P1(:,1:3)\a2h(u(1:2,:)); % direction vectors
                x2 = P2(:,1:3)\a2h(u(3:4,:));
                X = zeros(3,size(x1,2));
                for i=1:size(x1,2)
                    % n = cross(x2(:,i),x1(:,i)); slow
                    n = cr(x2(:,i),x1(:,i));
                    A = [x1(:,i) -x2(:,i) n];
                    if rank(A)>2
                        a = A\T;
                        X(:,i) = C1+a(1)*x1(:,i)+a(3)/2*n;
                    else
                        X(:,i) = inf(3,1);
                    end
                end
            otherwise
                error('unknown method');
        end
        % reprojections
        z = zeros(2,size(X,2));
        [v1,z(1,:)] = X2u(X,P1);
        [v2,z(2,:)] = X2u(X,P2);
        v = [v1(1:2,:);v2(1:2,:)];
        e = vnorm(v-u);
        % undo 3D normalization
        X = h2a(TP*a2h(X));
        a = PPX2ae(P1,P2,X);
        if ~isempty(a)
            a = a(3,:);
        end
    end
else
    X  = [-2  2  2 -2 -2  2  2 -2  
          -2 -2  2  2 -2 -2  2  2
           6  6  6  6  8  8  8  8];
    % simulate cameras  
    K{1} = [1000    0  500
             0   1000  500
             0      0    1];
    K{2} = K{1};
    R{1} = [1 0 0
            0 1 0
            0 0 1];
    R{2} = R{1};
    C{1} = [-1/2;0;0];
    C{2} = [ 1/2;0;0];
    P{1} = K{1}*[R{1} -R{1}*C{1}];
    P{2} = K{2}*[R{2} -R{2}*C{2}];
    % exact projections
    P = [P{1};P{2}];
    u = [h2a(X2u(X,P(1:3,:)));h2a(X2u(X,P(4:6,:)))];
    rng(0,'twister'); n = 2*(rand(size(u))-0.5);
    % large translation
    X0 = [1;pi*10^4;exp(1)*10^8];
    P1 = [K{1}*R{1}*[eye(3) -C{1}-X0];  K{2}*R{2}*[eye(3) -C{2}-X0]];
    X1 = X + X0*ones(1,size(X,2));
    [Xr,v]=uu2X(u,P1,'PROJ');
    test.fname = mfilename;
    test.test = {'Exact projections PROJ'};
    test.ok = max(vnorm(Xr-X1))<1e-8;
    %    
    [Xr,v]=uu2X(u,P1,'TRAN');
    test.test{end+1} = {'Exact projections TRAN'};
    test.ok(end+1) = max(vnorm(Xr-X1))<1e-8;   
    % large scale    
    P1 = [10^(-8)*[K{1}*R{1}*[eye(3) -C{1}]];  10^8*[K{2}*R{2}*[eye(3) -C{2}]]];
    X1 = X;
    [Xr,v]=uu2X(u,P1,'PROJ');
    test.test{end+1} = {'Large scale PROJ'};
    test.ok(end+1) = max(vnorm(Xr-X1))<1e-7;   
    %
    [Xr,v]=uu2X(u,P1,'TRAN');
    test.test{end+1} = {'Large scale TRAN'};
    test.ok(end+1) = max(vnorm(Xr-X1))<1e-7;   
    % test
    [Xr,v]=uu2X(u,P,'PROJ');
    test.test{end+1} = {'Test 1 PROJ'};
    test.ok(end+1) = max(vnorm(Xr-X1))<1e-10;   
    %
    [Xr,v]=uu2X(u,P,'TRAN');
    test.test{end+1} = {'Test 1 TRAN'};
    test.ok(end+1) = max(vnorm(Xr-X1))<1e-14;   
    %
    [Xr,v]=uu2X(u+n,P,'PROJ');
    test.test{end+1} = {'Test 2 PROJ'};
    test.ok(end+1) = max(vnorm(Xr-X))<0.07;   
    %    
    [Xr,v]=uu2X(u+n,P,'TRAN');
    test.test{end+1} = {'Test 2 TRAN'};
    test.ok(end+1) = max(vnorm(Xr-X))<0.07; 
    % multi view triangulation
    Pm = {K{1}*[eye(3) -[0;0;0]]; K{1}*[eye(3) -[5;0;0]]; K{1}*[eye(3) -[10;0;0]]; K{1}*[eye(3) -[50;0;0]]; K{1}*[eye(3) -[100;0;0]]};
    um = cellfun(@(p) h2a(X2u(X,p)),Pm,'UniformOutput',false);
    Pm = cat(1,Pm{:});
    um = cat(1,um{:});
    [Xr,v,z,e,a] = uu2X(um,Pm);
    test.test{end+1} = {'Multiview tringulation'};
    test.ok(end+1) = max(vnorm(Xr-X))<1e-15; 
    % error influence
    sX = max(vnorm(kron(X,ones(1,size(X,2)))-kron(ones(1,size(X,2)),X))); % max distance between points X
    M = cumprod([0.1 2*ones(1,11)]);
    e = zeros(2,size(M,2));
    a = PPX2ae(P(1:3,:),P(4:6,:),X); a = mean(a(3,:)); % mean appical angle at X
    for i=1:size(M,2)
        e(:,i) = [mean(vnorm(uu2X(u+M(i)*n,P,'PROJ')-X));mean(vnorm(uu2X(u+M(i)*n,P,'TRAN')-X))];
    end
    subfig(2,3,1); loglog(M*max(abs(n(:))),(e/sX),'.-');axis tight;
    xlabel('uniform image noise max');
    ylabel('3D error [% scene size]');
    legend('PROJ','TRAN');
    title(sprintf('3D reconstruction error, mean(apic\\angle(X)=%.1f deg',180*a/pi));
    % depth influence
    X = [zeros(2,80); 0.2:0.2:16];
    % lateral stereo
    u = [h2a(X2u(X,P(1:3,:)));h2a(X2u(X,P(4:6,:)))];
    u2 = u;
    u(3,:) = u(3,:)+1; % 1 pixel error along epipolar lines
    P2 = [P(1:3,:);KRC2P(K{2},a2r([0;1;0],2*atan(1/2000)),C{2})];
    % forward motion
    P3 = [KRC2P(K{1},R{1},[0;0;0]);KRC2P(K{2},R{2},[0;0;-1])];
    X3 = [0.5*ones(1,80);zeros(1,80);0.2:0.2:16];
    u3 = [h2a(X2u(X3,P3(1:3,:)));h2a(X2u(X3,P3(4:6,:)))];
    u3(3,:) = u3(3,:)+1; % 1 pixel error along epipolar lines
    e = [vnorm(uu2X(u,P,'PROJ')-X); vnorm(uu2X(u2,P2,'PROJ')-X); vnorm(uu2X(u3,P3,'PROJ')-X3)];
    a = [PPX2ae(P(1:3,:),P(4:6,:),X); PPX2ae(P3(1:3,:),P3(4:6,:),X3)];
    a = a([3 6],:); % appical angle at X
    subfig(2,3,2); plot(X(3,:),[100*e(1:2,:)./([1;1]*X(3,:));1.80*a(1,:)/pi],'.-');axis tight;
    xlabel('z coordinate');
    ylabel('3D rel err 100*e/z [%] & ap \angle/100 [\circ]');
    legend('100*e/z pix','100*e/z rot','ap \angle/100');
    title('Lateral Stereo (b=1m 1000px/60\circ)');
    subfig(2,3,3);
    plot(X3(3,:),100*e(3,:)./(X3(3,:)),'.-',X3(3,:),180*a(2,:)/pi,'.-r');axis tight;
    legend('e/z pix','ap \angle');
    xlabel('z coordinate');
    ylabel('3D rel err e/z & ap \angle [\circ]');
    legend('e/z pix','ap \angle');
    title('Forward Stereo (b=1m 1000px/60\circ)');
    % test output
    X = test;
end

function z = cr(x,y)
  z = [x(2)*y(3)-x(3)*y(2);-x(1)*y(3)+x(3)*y(1);x(1)*y(2)-x(2)*y(1)];
return
