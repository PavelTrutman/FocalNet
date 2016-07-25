% [o,A,x] = dtctspot(im,decim,dthrprc,minR,maxR,dThr,msk,doShow) - detect bright spots
%
% im     ... m x n x 1 image
% decim  ... 1:decim:end image decimation 
% dthrprc... thereshold 0<1  ... percentile 
%                       >=1 ... absolute value 
% minR   ... bright center outer radius
% maxR   ... dark perifery inner radius
% dThr   ... peak decision threshold
% msk    ... mask [[r;c], ...] for non-(max)ima supression
% doShow ... 1 = show images
% o      ... o{i} = [r;c;v;l] 
% A      ... A(i) = area
% x      ... x(:,i) = centroid

% (c) T.Pajdla, www.neovision.cz, Aug 4, 2005
function [o,A,x] = dtctspot(im,decim,dthrprc,minR,maxR,dThr,msk,doShow)

if dthrprc<1
    dthrprc = 100 * dthrprc;
    % decimate image to speed up
    imv = im(1:decim(1):end,1:decim(2):end);
    % the trace must brighter then the average intensity
    ix  = imv>mean(imv(:));
    % the trace is brighter than the dthrprc % of intensities
    thr = prctile(imv(ix),dthrprc);
else
    thr = dthrprc; 
end
% masks
% only bright pixels are interesting
imt = double(im>thr);
% ever used data: use separable filtering
v = conv2(imt,ones(8*maxR+3,1),'same');
v = conv2(v  ,ones(1,8*maxR+3),'same');
% zero parts that are never used for result to speed up convolutions
imv    = im;
v      = ~v;
imv(v) = 0;
% center: use separable filtering
ctrN = (2*minR+1)^2;
c  = conv2(imv,ones(2*minR+1,1),'same');
c  = conv2(c  ,ones(1,2*minR+1),'same');
c  = c/ctrN;
% background: use separable filtering
bckN = (4*maxR+1)^2;
b  = conv2(imv,ones(4*maxR+1,1),'same');
b  = conv2(b  ,ones(1,4*maxR+1),'same');
b  = b/bckN;
% hight center, low background
dcb = c - b;
imb = imt & dcb>dThr;
% label potential regions
[il,lN] = bwlabel(imb,8);
for i = 1:lN
    a(i) = length(il(il==i));
end
% centroids & areas
o = {}; A = []; x = [];
for j=1:lN
    [p,q] = find(il==j);
    ix = sub2ind(size(il),p,q);     
    iv = im(ix);
    o{j} = [p';q';iv';j*ones(1,size(p,1))];
    A(j) = length(p);
    x(:,j) = [mean(p);mean(q)];
end
% mask out those that are close to others
y = round(x);
if ~isempty(msk)
    for j=1:lN
        ix = sub2ind(size(il),msk(1,:)+y(1,j),msk(2,:)+y(2,j));
        rix(j) = ~any(imb(ix));
    end
    if ~isempty(o)
        o = o(rix);
        A = A(rix);
        x = x(:,rix);
    end
end
%% Plots
if doShow
    %% show image
    subfig(4,5,1);
    imagesc(im);
    axis on; colormap gray; axis image
    title('im - image')
    %% show thresholded image
    subfig(4,5,2);
    imagesc(imt);
    axis on; colormap gray; axis image
    title('imt - thresholded')
    %% show ever used pixels
    subfig(4,5,3);
    imagesc(v);
    axis on; colormap gray; axis image
    title('v - valid')
    %% center
    subfig(4,5,4);
    imagesc(c);
    axis on; colormap gray; axis image
    title('c - center')
    % background
    subfig(4,5,5);
    imagesc(b);
    axis on; colormap gray; axis image
    title('b - backround');
    % center - background
    subfig(4,5,6);    
    image(c-b);    
    axis on; colormap gray; axis image
    title('c-b');
    % center high & background low
    subfig(4,5,7);
    imagesc(imb);
    axis on; colormap gray; axis image
    title(sprintf('imb=(c-b)>%d & imt',dThr))
    % 
    subfig(4,5,8);
    imagesc(il);
    axis on; colormap gray; axis image
    title('il')
    %
    if lN>0
        imshow(label2rgb(il,@jet,'k'),'InitialMagnification','fit');        
    end
    title([num2str(lN) ' peak region(s)'])
    drawnow;
end
return