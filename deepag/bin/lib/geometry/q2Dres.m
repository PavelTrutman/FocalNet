% [r,y,w] = q2Dres(Q,x[,s,dbgPath]) - Conic residuals
%
% Q ... 3 x 3 conic matrix
% x ... 2 x n points
% s ... [n d] - iteration stops if it>n or no point moves more 
%               than d fraction of its distance ([10 0.01] implicitly)
%       0 < n < 1 ... for each x: angle(y-x,normal(y)) < n*360 degrees 
% dbg ... debug flag 
%         w{i} ... 2 y's for each iteration 
% r ... 1 x n residuals
% y ... the points on which r is realized
% 
% See q2Dfit, 
% Run q2Dres; for demo

% (c) T.Pajdla, www.neovision.cz, Sep 24, 2006
function  [r,x,w] = q2Dres(Q,x,s,dbgPath)
    w = [];
    if nargin<4, dbgPath=false; end % implicit no debug
    if nargin<3, s = [10 0.01] ; end 
    if nargin>0 % not demo
        if s(2)>0
            Q = Q/norm(Q(:),2); % normalization to alleviate numerical problems
            x = a2h(x); % homogeneous coordinates
            z = x; % keep the initial points
            ix = true(1,size(x,2)); % index of points to be updated
            if dbgPath
                k = 1;
                w{k} = x; k = k+1;
            end
            i = 1;
            while 1
                xi = x; % before the iteration
                %% the closest point to z on the intersection of the differential and the xy plane
                xx = x(:,ix); % to be updated
                e = sum(xx.*(Q*xx)); % values of the quadric at x
                J = 2*Q(1:2,:)*xx; % normal to the quadric at x projected to xy plane
                a = -e./sum(J.*J); % coefficient
                x(1:2,ix) = xx(1:2,:) + ([1;1]*a).*J;  % update x
                if dbgPath
                    w{k} = x; k = k+1;
                end
                %% the closest point to z on the intersection of the next differential at the current
                %% point
                xx = x(:,ix);
                J  = 2*Q(1:2,:)*xx; % normal to the quadric at x projected to xy plane
                J  = J./([1;1]*sqrt(sum(J.*J))); % normalize
                Jt = [0 -1;1 0]*J; % direction vector to the line pependicular to J at x, already normed
                x(1:2,ix) = Jt.*([1;1]*sum(Jt.*z(1:2,ix))) +  J.*([1;1]*sum(J.*xx(1:2,:))); % update x
                if dbgPath
                    w{k} = x; k = k+1;
                end
                %% stopping criterion
                if i>s(1), break; end
                if s(2)>0
                    dx = vnorm(x(1:2,ix)-xi(1:2,ix))./vnorm(x(1:2,ix)-z(1:2,ix));
                    ix(ix) = dx>s(2);
                    if ~any(ix), break; end
                end
                i = i + 1;
            end
        else
            Q = Q/norm(Q(:),2); % normalization to alleviate numerical problems
            x = a2h(x); % homogeneous coordinates
            z = x; % keep the initial points
            ix = true(1,size(x,2)); % index of points to be updated
            if dbgPath
                k = 1;
                w{k} = x; k = k+1;
            end
            for i=1:s(1)
                %% the closest point to z on the intersection of the differential and the xy plane
                e = sum(x.*(Q*x)); % values of the quadric at x
                J = 2*Q(1:2,:)*x; % normal to the quadric at x projected to xy plane
                a = -e./sum(J.*J); % coefficient
                x(1:2,:) = x(1:2,:) + ([1;1]*a).*J;  % update x
                if dbgPath
                    w{k} = x; k = k+1;
                end
                %% the closest point to z on the intersection of the next differential at the current
                %% point
                J  = 2*Q(1:2,:)*x; % normal to the quadric at x projected to xy plane
                J  = J./([1;1]*sqrt(sum(J.*J))); % normalize
                Jt = [0 -1;1 0]*J; % direction vector to the line pependicular to J at x, already normed
                x(1:2,:) = Jt.*([1;1]*sum(Jt.*z(1:2,:))) +  J.*([1;1]*sum(J.*x(1:2,:))); % update x
                if dbgPath
                    w{k} = x; k = k+1;
                end
            end
        end
        r = sqrt(sum((x-z).^2)); % the distance
    else % demo
        [U,D,V]=svd(rand(2));
        %V = eye(2);
        V=[V [0;0]; [0,0,1]]*[1 0 -3;0 1 -2;0 0 1];
        q=V'*diag([1 4*randn(1) -1])*V;
        x = 5*rand(2,100);
        %x=h2a(sampleq(q,100))+0.3;
        subfig(2,2,1);
        plot(x(1,:),x(2,:),'.'); hold on
        axis([-1 6 -1 6])
        %axis([-0 2.5 0 2.5])
        axis equal        
        tic; [d,y,w] = q2Dres(q,x,[10 0.01],true); toc
        plot(y(1,:),y(2,:),'.r');
        line([x(1,:);y(1,:)],[x(2,:);y(2,:)]);
        set(plotq(q),'linewidt',2,'color','g');
        subfig(2,2,2);
        plot(w{1}(1,:),w{1}(2,:),'.'); hold on
        for k=2:length(w)
            plot(w{k}(1,:),w{k}(2,:),'.'); hold on
            line([w{k-1}(1,:);w{k}(1,:)],[w{k-1}(2,:);w{k}(2,:)]);
        end
        set(plotq(q),'linewidt',2,'color','g');        
        axis([-1 6 -1 6])
        %axis([-0 2.5 0 2.5])
        axis equal
    end
return