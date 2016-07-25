% t = xconstr2T(x1,x2,d[,p]) - L_infty solution to camera translation for known rotation
%
% x1, x2 ... image points, 3 x n
% d      ... [d1 d2] = tolerance circle radii 
% sp     ... 0 (or mising) - d = radius in pixels in the 2D image  
%            1             - d = radius in radians on the sphere 
% p      ... true=plot intermediate results

% R. Hartley & F. Kahl. ICCV 2007.
%
% (c) T.Pajdla, pajdla@cmp.felk.cvut.cz
% 5 May 2008
function t = xconstr2T(x1,x2,dt,p,sp)

if nargin<5
    sp = 0;
end
if nargin<4
    p=false;
end
if length(dt)==1
    dt = [dt dt];
end

switch sp
    case 0 % in image 
        % normalize to the images
        x1 = a2h(h2a(x1));
        x2 = a2h(h2a(x2));
        w  = (x1+x2)/2;
        p1 = zeros(3,size(x1,2));
        p2 = zeros(3,size(x1,2));
        % plane translation constraints from images points
        for j=1:size(x1,2)
            n = [0 1 0;-1 0 0;0 0 0]*(x1(:,j)-w(:,j));
            z1 = x1(:,j)+dt(1)*n/sqrt(norm(n)^2-dt(1)^2);
            z2 = x1(:,j)-dt(1)*n/sqrt(norm(n)^2-dt(1)^2);
            p1(:,j) = cross(w(:,j),z1);
            p1(:,j) = p1(:,j)/norm(p1(:,j));
            p2(:,j) = cross(z2,w(:,j));
            p2(:,j) = p2(:,j)/norm(p2(:,j));
        end
        p1 = p1(:,all(imag(p1)==0));
        p2 = p2(:,all(imag(p2)==0));
    case 1 % on the sphere
        % normalize to the sphere
        x1 = x1./([1;1;1]*vnorm(x1));
        x2 = x2./([1;1;1]*vnorm(x2));
        x = cross(x1,x2);
        sa = vnorm(x); % sin alpha
        ix = sa>sin(dt(1)+dt(2)); % select active constraints
        sa = sa(ix);
        x = x(:,ix); 
        if sum(ix)>0
            x = x./([1;1;1]*sa);
            ca = sum(x1(:,ix).*x2(:,ix)); % cos aplha
            sb = sqrt(sin((dt(1))^2+2*sin(dt(1))*sin(dt(2))*ca+sin(dt(2))^2)./sa.^2); % sin beta
            cb = sqrt(1-sb.^2);
            w = (sin(dt(2))*x1(:,ix)+sin(dt(1))*x2(:,ix))./([1;1;1]*(sb.*sa));
            z = cross(x,w);
            p1 = ([1;1;1]*sb).*z+([1;1;1]*cb).*x;
            p1 = p1./([1;1;1]*vnorm(p1));
            p2 = ([1;1;1]*sb).*z-([1;1;1]*cb).*x;
            p2 = p2./([1;1;1]*vnorm(p2));
        else
            p1 = [];
            p2 = [];
        end
end
% solutions
pp = [p1 p2];
switch(size(pp,2))
    case 0
        t = inf;
    case 1
        t = pp;
    otherwise
    if size(pp,2)>100
        % find 4 feasible translations by LP
        pp = pp';
        b = zeros(size(pp,1),1);
        pp = [pp ;[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]];
        b = [b;[1 1 1 1 1 1]'];
        t(:,1) = linprog([0; 0; 1],-pp,b);
        t(:,2) = linprog([0; 0;-1],-pp,b);
        t(:,3) = linprog([0; 1; 0],-pp,b);
        t(:,4) = linprog([0;-1; 0],-pp,b);
    else % faster for small number of matches
        % all intersections
        % replicate pp's
        ix1 = reshape(ones(size(pp,2),1)*(1:size(pp,2)),1,[]);
        ix2 = reshape((1:size(pp,2))'*ones(1,size(pp,2)),1,[]);
        % n*(n-1)/2 pairwise intersections
        z = cross(pp(:,ix1),pp(:,ix2));
        z = z(:,vnorm(z)>eps);
        z = z./([1;1;1]*vnorm(z));
        if size(z,2)<10000 % faster but memory greedy
            % replicate intersections and pp's
            ix1 = reshape((1:size(pp,2))'*ones(1,size(z,2)),1,[]);
            ix2 = reshape(ones(size(pp,2),1)*(1:size(z,2)),1,[]);
            % all pairwise signs
            s = sum(pp(:,ix1).*z(:,ix2));
            s = reshape(s,size(pp,2),[]);
            ix = find(all(s>-10*eps)); % vertices of the convex hull of solutions
        else
            ix = false(1,size(z,2)); % all points
            for i=1:size(z,2)
                ix(i) = all((z(:,i)'*pp)>=-10*eps);
            end
            ix = find(ix);
        end
        t = z(:,ix);
    end
    t = t./(ones(size(t,1),1)*vnorm(t));
end
% plot
if p
    col = {'r';'y'};
    switch sp
        case 0
            dt = dt(1);
            % camera
            camplot(eye(3,4));
            if ~ishold, hold; end
            % sphere
            [X,Y,Z]=sphere(20);h=surf(X,Y,Z);
            fi = 0:0.1:2*pi;
            cx = [cos(fi);sin(fi);0*fi];
            set(h,'facecolor','g');set(h,'edgecolor','g'); grid on
            axis equal; xlabel('x'); ylabel('y'); zlabel('z');
            %
            for j=1:size(x1,2)
                % point constraints in images
                plot3(x1(1,j),x1(2,j),1,'.',x2(1,j),x2(2,j),1,'.');
                R = [1 0 x1(1,j);0 1 x1(2,j);0 0 1]; Q = inv(R)'*diag([1 1 -dt^2])*inv(R);
                y = a2h(h2a(sampleq(Q))); plot3(y(1,[1:end 1]),y(2,[1:end 1]),y(3,[1:end 1]),'b');
                R = [1 0 x2(1,j);0 1 x2(2,j);0 0 1]; Q = inv(R)'*diag([1 1 -dt^2])*inv(R);
                y = a2h(h2a(sampleq(Q))); plot3(y(1,[1:end 1]),y(2,[1:end 1]),y(3,[1:end 1]),'b');
            end
            for j=1:size(p1,2)
                view(2);
                [h,vx]=plotline(p1(:,j));if ~isempty(h),delete(h);plot3(vx(1,:),vx(2,:),[1 1]); end
                [h,vx]=plotline(p2(:,j));if ~isempty(h),delete(h);plot3(vx(1,:),vx(2,:),[1 1]); end
                view(3)
                % on the unit sphere
                X = 1.01*xy2R3([0;0;1],p1(:,j))*cx;
                plot3d(X(:,[1:end 1]),col{1});
                X = 1.01*xy2R3([0;0;1],p2(:,j))*cx;
                plot3d(X(:,[1:end 1]),col{2});
            end
        case 1
            if ~ishold, hold on; end
            % sphere
            [X,Y,Z]=sphere(20);h=surf(X,Y,Z);
            fi = 0:0.1:2*pi;
            cx = [cos(fi);sin(fi);0*ones(size(fi))];
            set(h,'facecolor','g');set(h,'edgecolor','g'); grid on
            axis equal; xlabel('x'); ylabel('y'); zlabel('z');
            %
            y1 = sampleq(diag([1 1 -dt(1)^2]));
            y2 = sampleq(diag([1 1 -dt(2)^2]));
            for j=1:size(x1,2)
                % point constraints on the sphere
                R = xy2R3([0;0;1],x1(:,j));
                y = 1.01*R*y1;
                plot3d(y(:,[1:end 1]),'b');
                R = xy2R3([0;0;1],p1(:,j));                
                X = 1.01*R*cx;
                plot3d(X(:,[1:end 1]),col{1});                
                R = xy2R3([0;0;1],x2(:,j));
                y = 1.01*R*y2;
                plot3d(y(:,[1:end 1]),'b');
                R = xy2R3([0;0;1],p2(:,j));                                
                X = 1.01*R*cx;
                plot3d(X(:,[1:end 1]),col{2});                
            end           
    end
    plot3d(p1,'r.');
    plot3d(-p2,'y.');
    if size(t,1)==3
        for i=1:size(t,2),
            plot3([0 1.3*t(1,i)],[0 1.3*t(2,i)],[0 1.3*t(3,i)],'linewidth',2);
        end
        view([120 30]);
        % angles
        ix1 = reshape((1:size(t,2))'*ones(1,size(t,2)),1,[]);
        ix2 = reshape(ones(size(t,2),1)*(1:size(t,2)),1,[]);
        a = max(acos(sum(t(:,ix1).*t(:,ix2))));
        title(num2str(180/pi*a));
    end
end