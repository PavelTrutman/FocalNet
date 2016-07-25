% [H,d] = xy2H(x,y[,m,dp,op]) - Homography estimation
% 
% \exists a \in R such that a y ~ H x where ~ stands for "approximates"
%
% x ... 3 x n 2D points 
% y ... 3 x n 2D points
% H ... 3 x 3 homography matrix
% d ... 1 x n residuals computed by using resid
% m ... {method, resids}
%        method = 'DLT' ... direct liear normalized transform (implicit)
%                 'GST' ... DLT followed by the the gold standard 
%                           ~ full non-linear bundle adjustment of sum(d.^2)
%        resids = 'x-x' ... error in x
%                           dx = sqrt(sum((h2a(inv(H)*y)-h2a(x)).^2));
%                 'y-y' ... error in y
%                           dy = sqrt(sum((h2a(y)-h2a(H*x)).^2)); 
%                 'x-y' ... mean odf the error in x + error y (implicit)
%                           d  = (dx + dy)/2;
%         
% dp ... debug info & plots flag = 0 ... do nothing (implicit)
%                                  1 ... print & plot
%               
% op ... optimization parameters
%        m = 'DLT' ... ignored
%            'GST' ... opts = optimset(...)
%                      implicit: optimset('Display','none','TolX',1e-5,'MaxIter',100);
%
% R.Hartley, A.Zisserman. Multiple View Geometry. 2nd edition. Cambridge
% Press 2003.

% (c) T.Pajdla, www.neovision.cz, Nov 1 2005
function [H,d] = xy2H(x,y,mth,dp,opts)
%
if nargin==0
    xy2Hdemo;
    return
end
%
if nargin<5
    opts  = optimset('Display','none','TolX',1e-5,'MaxIter',100);
end
%
if nargin<4
    dp = 0;
end
if nargin<3
    mth = {'DLT','x-y'};
end
% upgrade to homogeneous coordinates if necessary
if size(x,1)<3
    x = a2h(x);
end
if size(y,1)<3
    y = a2h(y);
end
%
if strcmp(mth{1},'DLT') || strcmp(mth{1},'GST')
    % Tx = x2nx(x);       % normalization
    % Ty = x2nx(y);
    [Tx,Ty] = xy2nxy(x,y);  % simultaneous normalization makes a better sense for the same cameras
    xn = Tx*x;            % normalized points
    yn = Ty*y;        
    %   A  = zeros(3*size(x,2),9);  % construct the design matrix        
    %   for i=1:size(x,2)
    %         % Intuitive
    %         % A([1:3]+3*(i-1),:) = kron([ 0       -yn(3,i) yn(2,i);...
    %         %                             yn(3,i)  0      -yn(1,i);...
    %         %                            -yn(2,i)  yn(1,i) 0      ],xn(:,i)');
    %         % Efficient
    %         A([1:3]+3*(i-1),:) = ...
    %              [[0;yn(3,i);-yn(2,i)]*xn(:,i)',...
    %               [-yn(3,i);0;yn(1,i)]*xn(:,i)',...
    %               [yn(2,i);-yn(1,i);0]*xn(:,i)'];
    %     end
    % Even more efficient
    A = [[ zeros(size(yn,2),3) -yn([3 3 3],:)' yn([2 2 2],:)'].*xn([1 2 3 1 2 3 1 2 3],:)'
         [ yn([3 3 3],:)' zeros(size(yn,2),3) -yn([1 1 1],:)'].*xn([1 2 3 1 2 3 1 2 3],:)'
         [-yn([2 2 2],:)' yn([1 1 1],:)' zeros(size(yn,2),3) ].*xn([1 2 3 1 2 3 1 2 3],:)'];         
    % balancing has mized effect since it may blow up the small values
    % A = A./((vnorm(A')'+1)*ones(1,size(A,2)));
    % use eig instead of [U,S,V] = svd(A,'econ');
    if 0
        [U,S,V] = svd(A,'econ');
        H = reshape(V(:,end),3,3)';
    else
        [V,D] = eig(A'*A); % faster but maybe less robust
        H = reshape(V(:,1),3,3)';
    end
    H = inv(Ty)*H*Tx;   % denormalization
    H = H/norm(H);
    d = xy2Hres(H,x,y,mth{2});
end
if strcmp(mth{1},'GST')
    if dp % set debug info
        opts.Display = 'iter'; 
    end
    h = reshape(H,9,1); % to a vector
    h = h/norm(h);  % unit length
    %h = x2sphp(h);  % spherical parametrization by 8 parameters
    h = lsqnonlin(@xy2Hres,h,[],[],opts,x,y,mth{2}); % optimize
    %h = sphp2x(h);  % unit vector
    h = h/sqrt(sum(h.^2));
    d = xy2Hres(h,x,y,mth{2});        
    H = reshape(h,3,3);
end
% debug plots
if dp > 0
    f(1) = subfig(3,5,1);
    plot(x(1,:),x(2,:),'b.');hold;
    plot(y(1,:),y(2,:),'r.');
    title('points: b. = x, r. = y');
    %
    f(2) = subfig(3,5,2);hold;
    plot(xn(1,:),xn(2,:),'b.');
    title('norm. points: b. = x');
    %
    f(3) = subfig(3,5,3);hold;
    plot(yn(1,:),yn(2,:),'r.');
    title('norm. points: r. = y');
    %
    f(4) = subfig(3,5,4);hold;
    xy = h2a(inv(H)*y);
    dx = xy-h2a(x);
    quiver(x(1,:),x(2,:),dx(1,:),dx(2,:),0);
    plot(x(1,:),x(2,:),'b.');
    plot(xy(1,:),xy(2,:),'c.');
    title('mapped points: b. = x');
    %
    f(5) = subfig(3,5,5);hold;
    yx = h2a(H*x);
    dy = yx-h2a(y);
    quiver(y(1,:),y(2,:),dy(1,:),dy(2,:),0);
    plot(y(1,:),y(2,:),'r.');    
    plot(yx(1,:),yx(2,:),'m.');
    title('mapped points: r. = y');
    %
    f(6) = subfig(3,5,6);hold;
    [hn,hx]=hist(dx',25);
    stairs([hx hx],hn)    
    title('hist dx')
    %
    f(7) = subfig(3,5,7);hold;
    [hn,hy]=hist(dy',25);
    stairs([hy hy],hn)    
    title('hist dy')
    drawnow
end
return
%% homography residuals in the spherical parametrization
function e = xy2Hsres(p,x,y,obj)
    H = sphp2x(p);
    e = xy2Hres(H,x,y,obj);
return
%% demo
function xy2Hdemo
    opts  = optimset('Display','none','TolX',1e-5,'MaxIter',20);
    %rand('state',313971);
    %randn('state',313971);    
    % points & H
    H = [1.0            0.0*randn(1)   0.0*randn(1); 
         0.0*randn(1)   1.0            0.0*randn(1);
         0.05*(rand(1))   0.0*randn(1)              1];
    [x,y] = meshgrid(1:10:100,1:10:100);
    x = [x(:)';y(:)'];
    y = h2a(H*a2h(x));
    % noise
    xx= x + 0.3*randn(size(x));
    yy= y + 0.3*randn(size(y));
    % get H by the DLT fit and plot debug results
    [H1,d1]=xy2H(xx,yy,{'DLT','y-y'},1);
    % get H by the gold standard fit
    [H2,d2]=xy2H(xx,yy,{'GST','y-y'},1,opts);   
    %
    subfig(3,2,4);
    plot(1:length(d1),d1,'b.-',1:length(d2),d2,'r.-');
    title(sprintf('|residuals|: b - DLT, r - GST mean([DLT GST]) = [%.3f %.3f]',mean(d1),mean(d2)));
    subfig(3,2,6);
    plottab([H/H(3,3),H1/H1(3,3)-H/H(3,3),H2/H2(3,3)-H/H(3,3)],'',['  ';'H ';'  ';'H ';'- ';'H1';'H ';'- ';'H2']);
    % measure efficiency
    drawnow
    tic;
    for i=1:100
        [H1,d1]=xy2H(xx,yy,{'DLT','y-y'});
    end
    disp(sprintf('xy2H(xx,yy,''DLT'') took %.4f [s]',toc/100));
    tic;
    for i=1:1
        [H1,d1]=xy2H(xx,yy,{'GST','y-y'});
    end
    disp(sprintf('xy2H(xx,yy,''GST'') took %.4f [s]',toc/1));    
return
