% [io,A,ri,ro,roi,xi,xo,xoi] = imHim(im,H[,rs,rc]) - Homography mapping of an image
%
% H  ... 3 x 3 homography matrix w.r.t. coordinates [col row 1]
% im ... input image
% rs ... resolution = number of pixels for the smallest size of the output bounding box
% rc ... [cmin cmax rmin rmax] rectangle (full image implicitly)
% io ... output image
% A  ... A*H maps input to output image, i.e. inv(A*H) maps detections 
%        from the output image into th einput image  
% ri ... corners of input rectangle
% ro ... corners of output rectangle, h2a(H*a2h(ri)) 
% roi... corners of output rectangle in output image coordinates, h2a(A*a2h(ro))
% xi ... input image sample points
% xo ... ouput image sample points, xo = h2a(H*a2h(xi))
% xoi... output image sample point indices, i.e. xo in output image coordinates
%
% Pseudo-algorithm:
% [u;v] = bounding_box(h2a(H*rc_corners));
% [x;y] = h2a(inv(H)*a2h([u;v]));
% io(1:rnum,1:cnum) = interp2(im,x,y,'cubic'); 

% T.Pajdla, pajdla@cvut.cz, 2016-07-27
function [io,A,ri,ro,roi,xi,xo,xoi] = imHim(im,H,rs,rc)
if nargin>0
    if nargin<4 || isempty(rc)
        rc = size(im); rc = [1 rc(2) 1 rc(1)];
    end
    if nargin<3 || isempty(rs)
        rs = min([rc(2)-rc(1)+1 rc(4)-rc(3)+1]); 
    end
    % corners of the input region
    ri = [rc(1) rc(2) rc(2) rc(1)
          rc(3) rc(3) rc(4) rc(4)];
    % corners of the output region
    ro = h2a(H*a2h(ri));
    % bounding box of the output region
    bo = [min(ro(1,:)) max(ro(1,:)) min(ro(2,:)) max(ro(2,:))];
    % output image size
    sz = rs/min((bo(2)-bo(1)+1),(bo(4)-bo(3)+1));
    sz = [bo(2)-bo(1)+1 bo(4)-bo(3)+1]*sz; 
    % transform to output image pixel coordinates
    xbo = [bo(1) bo(2) bo(2) bo(1)
           bo(3) bo(3) bo(4) bo(4)]; 
    xboi = [1 sz(1) sz(1) 1  
            1 1     sz(2) sz(2)];
    A = xy2H(a2h(xbo),a2h(xboi));
    roi = h2a(A*a2h(ro));
    % output points
    % output image point ranges
    xo = {bo(1):((bo(2)-bo(1))/(sz(1)-1)):bo(2)
          bo(3):((bo(4)-bo(3))/(sz(2)-1)):bo(4)};
    [x,y] = meshgrid(xo{1},xo{2});
    xo  = [x(:)';y(:)']; 
    szo = [size(x,2) size(x,1)]; % output image size
    [xoi,yoi] =  meshgrid(1:szo(1),1:szo(2));
    xoi = [xoi(:)';yoi(:)']; clear yoi x y
    % input points
    xi = h2a(H\a2h(xo));
    % remove points outside the input rectangle
    ix = xi(1,:)>=rc(1) & xi(1,:)<=rc(2) & xi(2,:)>=rc(3) & xi(2,:)<=rc(4);
    xi = xi(:,ix);
    xo = xo(:,ix);
    xoi = xoi(:,ix);
    % create output image
    io0 = nan(szo(2),szo(1));
    ixo = sub2ind(size(io0),xoi(2,:),xoi(1,:));
    for i=1:size(im,3)
        ioi = io0; 
        ioi(ixo) = interp2(im(:,:,i),xi(1,:),xi(2,:),'linear'); 
        io(:,:,i) = ioi; 
    end
else % unit tests
    im = checkerboard(3,3,4);
    [c,r] = meshgrid(1:size(im,2),1:size(im,1));
    im = 10*im+r; 
    % Test 1 - identity
    H = eye(3);
    [o,A,ri,ro,roi,xi,xo,xoi] = imHim(im,H);
    subfig(3,4,1); imagescax(im); hold on; plot3d(ri(:,[1 2 3 4 1]),'.-r'); plot3d(xi,'.y'); title('1: input');
    subfig(3,4,2); imagescax(o); hold on; plot3d(roi(:,[1 2 3 4 1]),'.-r');  plot3d(xoi,'.y'); title('1: output');
    if size(xi)==size(xo)
        io(1) =  all(vnorm(xi-xo)<eps);
    else
        io(1) = false;
    end
    % Test 2 - intermal computation check
    H = [1 0 0;0 1 0;1/20 1/10 1];
    [o,A,ri,ro,roi,xi,xo,xoi] = imHim(im,H,100);
    subfig(3,4,1); imagescax(im); hold on; plot3d(ri(:,[1 2 3 4 1]),'.-r'); title('1: input');
    subfig(3,4,2); imagescax(o); hold on; plot3d(roi(:,[1 2 3 4 1]),'.-r'); title('1: output');
    if size(xi)==size(xo)
        io(2) =  all(vnorm(h2a(A*a2h(xo))-xoi)<1e-12);
    else
        io(2) = false;
    end
    % Test 3 - almost affin mapping
    [c,r] = meshgrid(1:360,1:240);
    im = 200*single(((c-180).^2+(r-130).^2-50^2)<=0);
    im = imfilter(im,fspecial('gaussian',21,3),'same');
    H = [0.99 1 10.1;0.1 2 0.65;1/1000 0 1];
    rc = [100 250 50 200];
    [o,A,ri,ro,roi,xi,xo,xoi] = imHim(im,H,2000,rc);
    spx = [0.5 3 30 2 1]; % ellipse subpixel detection parameters
    [p{1},q{1},~,~,~,f] = detectellipse(uint8(im),[],[],[],[],[],[],[],[],spx,true);
    [p{2},q{2},e,a,b] = detectellipse(uint8(o),[],6,0.01*numel(o),0.5*numel(o),[],[],[],[],spx,true);
    G = A*H;
    S = G'*q{2}*G;
    [s,B] = q2par(S);     
    x = plotq(S,[]);
    figure(f{2}); plot3d(x,'m.-');plot3d(s(1:2),'.m');title(sprintf('center error = %.3f',vnorm(p{1}(1:2)-s(1:2))));
    if size(xi)==size(xo)
        io(3) =  vnorm(p{1}(1:2)-s(1:2))<0.01;
    else
        io(3) = false;
    end
    % Test 4 - small shifts
    clear s p q e
    [c,r] = meshgrid(1:360,1:240);
    im = 100*single(((c-180).^2+(r-130).^2-20^2)<=0);
    im = imfilter(im,fspecial('gaussian',21,3),'same');
    spx = [2 5 30 6 4]; % ellipse subpixel detection parameters
    [p{1},q{1},~,~,~,f] = detectellipse(uint8(im),[],[],[],[],[],[],[],[],spx,true);
    for i=1:3
        H = [1 0 3*i/7;0 1 0;0 0 1];
        rc = [100 250 50 200];
        [o,A,ri,ro,roi,xi,xo,xoi] = imHim(im,H,2000,rc);
        [p{1+i},q{1+i}] = detectellipse(uint8(o),[],6,0.01*numel(o),0.5*numel(o),[],[],[],[],spx,true);
        G = A*H;S = G'*q{1+i}*G;[s{i},B] = q2par(S);
        e(i) = vnorm(p{1}(1:2)-s{i}(1:2));
    end
    io(4) =  all(e<0.01);
    % Test 5 - perspeective mappintgs
    clear s p q e
    [c,r] = meshgrid(1:360,1:240);
    im = 100*single(((c-180).^2+(r-130).^2-50^2)<=0);
    im = imfilter(im,fspecial('gaussian',21,2),'same');
    [p{1},q{1},~,~,~,f] = detectellipse(uint8(im),[],[],[],[],[],[],[],[],spx,4);
    e = [];
    if ~isempty(q{1})
        e(1) = 0;
    end
    for i=1:3
        H = [1 0 0;0 1 0;i/1000 i/1000 1];
        rc = [120 230 70 180];
        [o,A,ri,ro,roi,xi,xo,xoi] = imHim(im,H,2000,rc);
        [p{1+i},q{1+i}] = detectellipse(uint8(o),[],12,0.01*numel(o),0.5*numel(o),[],[],[],[],spx,4);
        if ~isempty(q{1+i})
            G = A*H;S = G'*q{1+i}*G;[s{i},B] = q2par(S);
            e(i+1) = vnorm(p{1}(1:2)-s{i}(1:2));
        end
    end
    subfig(3,2,2);
    ix = ~cell2mat(cellfunuo(@isempty,p));
    p = p(ix);
    e = e(ix);
    if ~isempty(p)
        ex = sort(cell2mat(cellfunuo(@(c) c(3:4),p(2:end))));
        plot(max(ex)./min(ex),e(2:end),'.-'); xlabel('excetricity');ylabel('error');title('perspectivity');
    end
    if isempty(e)
        io(5) = false;
    else
        io(5) =  all(e<0.02);
    end
end
return
