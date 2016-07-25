% [P,e,T] = uX2P(u,X[,cType,E,T]) - camera projection matrix P from 3D points X and their projections u
%
% cType = 'perspective' (implicit) or 'affine' ... estimates the standard cameras
%
%            u     ... 3 x n image points, affine coords for finite points
%            X     ... 4 x n homogeneous coordinates of 3D point X
%            P     ... 3 x 4 matrix
%            e     ... 2 x n residuals
%
%         'affine-E' ... estimates P & T from known u{i}, E{i}, and X{i}
%
%            u{i} = P * E{i} * T * X{i}
%
%            u{i}  ... 3 x n image points, affine coords for finite points
%            X{i}  ... 4 x n homogeneous coordinates of 3D point X
%            E{i}  ... 4 x 4 i-th known Euclidean transform
%            P     ... 3 x 4 matrix
%            e     ... 2 x n residuals
%            T     ... 4 x 4 Euclidean transoform

% (c) T.Pajdla, www.neovision.cz, Dec 5, 2004
function [P,e,T] = uX2P(u,X,ct,E,S)
if nargin>0
    if nargin<3
        ct = 'perspective';
    end
    switch ct
        case 'perspective'
            A           = eye(4);
            A(:,4)      = [-mean(X(1:3,:)')';1];
            Y           = A*X;
            s           = 1/mean(mean(sqrt(Y(:,1:3).^2)));
            A           = diag([s s s 1])*A;
            % A            = eye(4);
            Y           = A*X;            
            B           = eye(3);
            B(:,3)      = [-mean(u(1:2,:)')';1];
            v           = B*u;
            t           = 1/mean(mean(sqrt(v(:,1:2).^2)));
            B           = diag([t t 1])*B;
            v           = B*u;
            C = [];
            for i=1:size(v,2)
                z = xx(v(:,i));
                C = [C; [Y(1,i)*z Y(2,i)*z Y(3,i)*z Y(4,i)*z]];
            end
            [U,D,V] = svd(C,'econ');
            p = V(:,end);
            [pm,pmi]= max(abs(p));
            p = p/p(pmi);
            Q = reshape(p,3,4);
            P = B\(Q*A);
            P = P/vnorm(P(3,1:3));
            e = h2a(P*X)-h2a(u);
            T = eye(4);
        case 'affine'
            A           = eye(4);
            A(:,4)      = [-mean(X(1:3,:)')';1];
            Y           = A*X;
            s           = 1/mean(mean(sqrt(Y(:,1:3).^2)));
            A           = diag([s s s 1])*A;
            Y           = A*X;
            p1 = pinv(Y')*u(1,:)';
            p2 = pinv(Y')*u(2,:)';
            P  = [p1' ; p2' ; 0 0 0 1];
            P  = P*A;
            P  = P/vnorm(P(3,1:3));
            e  = h2a(P*X,zeros(1,size(X,2)))-h2a(u);
            T  = eye(4);
        case 'affine-E'
            % the normalization transformation of X
            A           = eye(4);
            Z           = cat(2,X{:});
            A(:,4)      = [-mean(Z(1:3,:)')';1];
            Y           = A*Z;
            s           = 1/mean(mean(sqrt(Y(:,1:3).^2)));
            A           = diag([s s s 1])*A;
            clear Y Z s
            % normalize
            for i=1:length(X)
                Y{i}    = A*X{i};
                n{i}    = u{i}(1:2,:)*pinv(Y{i});
            end
            B = []; C = [];
            for i=1:length(n)
                b = [kron(eye(4),n{i}(1,1:3)) -E{i}(1:3,:)' -[0;0;0;1]    zeros(4,4)
                    kron(eye(4),n{i}(2,1:3)) zeros(4,4)   -E{i}(1:3,:)' -[0;0;0;1]];
                c = [0;0;0;-n{i}(1,4); 0;0;0;-n{i}(2,4)];
                B = [B; b];
                C = [C; c];
            end
            s = pinv(B)*C;
            t = [reshape(s(1:12),3,4); 0 0 0 1];
            p = reshape(s(13:end),4,2)';
            % denormalize
            T = inv(t)*A;
            P = [p; 0 0 0 1];
            % evaluate the residual errors
            for i=1:length(X)
                v{i} = h2a(P*E{i}*T*X{i},zeros(1,size(X{i},2))); % projected points
                e{i} = v{i}-u{i}(1:2,:); % residuals
            end
            if 0 % plot reprojections
                uu = cat(2,u{:});
                vv = cat(2,v{:});
                plot(uu(1,:),uu(2,:),'r.',vv(1,:),vv(2,:),'g.');
                title(sprintf('error = %.3f',sum(sqrt(sum((uu-vv).^2)))/size(uu,2)))
                pause
            end
    end
else
    % 
    test.fname = mfilename;
    test.test = {'Demo'};
    test.ok = true;
    % an example
    s0      = 1; % 0 = no perturbation
    RndSeed = 137;
    
    X0   = a2h(trandn([3,36],2,-2,2)); % 3D points
    T    = tabc2tr([-2 -1 0 1 2 3]); % Camera position
    P    = [42 0  0  20 % Projection matrix
        0  42 10 0];
    for i=1:4
        E{i} = tabc2tr([rand(1,3) 60*rand(1,3)]);  % view positions
        X{i} = X0; % same points for all views
        u0{i}= P*E{i}*T*X{i}; % noise-free image points
        randn('state',RndSeed); % same pseudorandom noise
        u{i} = u0{i}+trandn([2 size(X{i},2)],1,-2,2); % noisy image points
    end
    XX  = cat(2,X{:}); % plot all 3D points
    f1 = subfig(3,4,1);
    plot3(XX(1,:),XX(2,:),XX(3,:),'.'); axis('equal');
    uu0 = cat(2,u0{:}); % all noise-free image points
    uu1 = cat(2,u{:}); % all noisy image points
    ct  = 'MonteCarlo';
    
    for i=1:length(u)
        u{i} = u{i}(1:2,:);
    end
    X0  = cat(2,X{:});
    XX  = X0;
    u0  = u;
    uu0 = cat(2,u{:});
    uu1 = uu0;
    P   = zeros(2,4);
    T   = eye(4);
    
    mu{1} = u; % the first estimate of image points = noisy image points
    for m=2:2
        uu = []; ee = []; pp = []; tt = []; ww = [];
        for k=1:400 % Monte Carlo Smapling
            for i=1:length(u)
                v{k}{i} = mu{m-1}{i}+s0*trandn([2 size(X{i},2)],2,-20,20); % perturb the current estimate of the image points
            end
            [p{k},e{k},t{k}] = uX2P(v{k},X,'affine-E',E); % estimate the best possible affine camera + the view positions
            for i=1:length(v{k})
                w{k}{i} = v{k}{i}+e{k}{i}; % projections of noisy points by the estimated camera from the estimated positions
            end
            for i=1:length(E)
                ep{k}{i}= p{k}(1:2,:)*E{i}*t{k}*X{i}-u0{i}; % differences of the projections of the correct points by the estimated camera from the correct image projections
            end
            uu   = [uu cat(2,v{k}{:})]; % all perturbed image projections
            ww   = [ww cat(2,w{k}{:})]; % all estimated image projections
            ee   = [ee cat(2,e{k}{:})]; % all residuals
            p1   = p{k}(1:2,:);
            t1   = t{k}(1:3,:);
            pp   = [pp ; p1(:)']; % all camera matrices
            tt   = [tt ; t1(:)']; % all view positions
            disp(num2str(k));
        end
        for i=1:length(w)
            W(:,:,i) = cat(2,w{i}{:}); % all sampled image projections
        end
        for i=1:length(e)
            me(i) = mean(sqrt(sum(cat(2,e{i}{:}).^2))); % mean squared error of the camera + the positions fit
        end
        mw   = mean(W,3); % the mean sampled image projections
        mpp  = mean(pp); % the mean camera matrices
        mtt  = mean(tt); % the mean view poitions
        mu{m}= mat2cell(mw,size(w{1}{1},1),cellfun('size',w{1},2)); % store estimated mean image projections
    end
    for i=1:length(mu)
        mue(i) = mean(sqrt(sum((cat(2,mu{i}{:})-cat(2,u0{:})).^2))); % mean squared error w.r.t. the correct image points
    end
    for i=1:length(ep)
        mep(i) = mean(sqrt(sum(cat(2,ep{i}{:})).^2)); % mean image projection errors for all samples
    end
    f2 = subfig(3,4,2);
    plot(uu(1,:),uu(2,:),'.'); axis('equal'); hold % perturbed image projections
    plot(ww(1,:),ww(2,:),'r.'); % estimated image projections
    plot(uu0(1,:),uu0(2,:),'g.','markersize',20); % noise-free image points
    plot(uu1(1,:),uu1(2,:),'c.','markersize',20); % noisy image points
    plot(mw(1,:),mw(2,:),'m.','markersize',20); % mean sampled image projections
    title('g - true; b,c - smpls; r,m - est');
    
    [sev,sei] = sort(me);
    if 0
        for i=1:size(W,2)
            plot(squeeze(W(1,i,sei)),squeeze(W(2,i,sei)),'r-');
            plot(squeeze(W(1,i,sei(1))),squeeze(W(2,i,sei(1))),'r.','markersize',20);
        end
    end
    f3 = subfig(3,4,3);
    [n,x] = hist(ee'); stairs(x,n); % residuals
    title('cam+pos fit residuals')
    f4 = subfig(3,4,4);
    plottab([reshape(mpp,2,4)-P;reshape(mtt,3,4)-T(1:3,:)],['dP';'  ';'dT';'  ';'  ';'  ']); % estimated R, T difference
    title('estimated R & T difference')
    f5 = subfig(3,4,5); % estimated & true camera matrices
    for i=1:size(pp,2)
        subplot(size(pp,2),1,i);
        [n,x] = hist(pp(:,i),10);  stairs(x,n,'r'); hold('on');
        plot([1 1]*mpp(i),[0 30],'m','linewidth',3);
        plot([1 1]*P(i),[0 30],'g','linewidth',3);
        axis([-inf inf 0 30]);
    end
    title('P')
    f5 = subfig(3,4,6); % estimated & true view positions
    t1 = T(1:3,:);
    for i=1:size(tt,2)
        subplot(size(tt,2),1,i);
        [n,x] = hist(tt(:,i),10);  stairs(x,n,'r'); hold('on');
        plot([1 1]*mtt(i),[0 30],'m','linewidth',3);
        plot([1 1]*t1(i),[0 30],'g','linewidth',3);
        axis([-inf inf 0 30]);
    end
    title('T');
    f6 = subfig(3,4,7)
    plot(mep);
    title('smpl mean errors');
    P = test;
end

