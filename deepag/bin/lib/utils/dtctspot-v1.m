% [iml,l,imt] = dtctspot(im,decim,dthrprc,minR,maxR,dThr,doShow) - detect bright spots
%
% im     ... m x n x 1 image
% decim  ... 1:decim:end image decimation 
% dthrprc... thereshold 0<1  ... percentile 
%                       >=1 ... absolute value 
% minR   ... bright center outer radius
% maxR   ... dark perifery inner radius
% dThr   ... peak decision threshold
% doShow ... 1 = show images
% iml    ... labeled image 
% l      ... labels \ {0}
% imt    ... thresholded image

% (c) T.Pajdla, www.neovision.cz, Aug 4, 2005
function [iml,l,imt] = dtctspot(im,decim,dthrprc,minR,maxR,dThr,doShow)

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

% only bright pixels are interesting
imt = double(im>thr);
% masks
ctrN = (2*minR+1)^2;
bckN = (4*maxR+1)^2;
% prfN = (4*maxR+1)^2 - (2*maxR+1)^2 = 
%      = 4^2*maxR^2 + 4*2*maxR + 1 - 2^2*maxR^2 - 2*2*maxR - 1
%      = (4^2-2^2)*maxR^2 + (4*2-2*2)*maxR 
%      = 3*4*maxR^2 + 2*2*maxR 
prfN = 4*3*maxR^2 - 4*maxR;
% ever used data: use separable filtering
% v   = conv2(imt,vldm,'same');
v = conv2(imt,ones(8*maxR+3,1),'same');
v = conv2(v  ,ones(1,8*maxR+3),'same');
% zero parts that are never used for result to speed up convolutions
imv    = im;
v      = ~v;
imv(v) = 0;
% center: use separable filtering
% c  = conv2(imv,ctrm,'same')/sum(ctrm(:));
c  = conv2(imv,ones(2*minR+1,1),'same');
c  = conv2(c  ,ones(1,2*minR+1),'same');
c  = c/ctrN;
% background: use separable filtering
% b  = conv2(imv,bckm,'same')/sum(bckm(:));
b  = conv2(imv,ones(4*maxR+1,1),'same');
b  = conv2(b  ,ones(1,4*maxR+1),'same');
b  = b/bckN;
% background: use separable filtering
% b  = conv2(imv,bckm,'same')/sum(bckm(:));
%mp = conv2(imv,ones(2*maxR+1,1),'same');
%mp = conv2(mp, ones(1,2*maxR+1),'same');
%p  = (bckN*b - mp)/prfN;
% hight center, low background
dcb = c - b;
imb = dcb>dThr;
cb  = double(imb);
% background: use separable filtering
% b  = conv2(imv,bckm,'same')/sum(bckm(:));
bcb  = conv2(cb,ones(4*maxR+1,1),'same');
bcb  = conv2(bcb,ones(1,4*maxR+1),'same');
% perifery: use separable filtering
% p  = conv2(bcb,prfm,'same')/sum(prfm(:));
mcb  = conv2(cb,ones(2*maxR+1,1),'same');
mcb  = conv2(mcb,ones(1,2*maxR+1),'same');
pcb  = bcb - mcb;
% potential peaks
pcb1      = pcb>0; 
imb(pcb1) = 0;
%imb = imb & imt; 
% label potential regions
iml = bwlabel(imb,8);
% labels \ {0} 
l   = unique(iml(iml(:)~=0));

if doShow
    %% show image
    subfig(4,5,1);
    imagesc(im);
    axis off; colormap gray; axis image
    title('im - image')
    %% show thresholded image
    subfig(4,5,2);
    imagesc(imt);
    axis off; colormap gray; axis image
    title('imt - thresholded')
    %% show ever used pixels
    subfig(4,5,3);
    imagesc(v);
    axis off; colormap gray; axis image
    title('v - valid')
    %% center
    subfig(4,5,4);
    imagesc(c);
    axis off; colormap gray; axis image
    title('c - center')
    % background
    subfig(4,5,5);
    imagesc(b);
    axis off; colormap gray; axis image
    title('b - backround');
    % center - background
    subfig(4,5,6);    
    image(c-b);    
    axis off; colormap gray; axis image
    title('c-b');
    % center high & background low
    subfig(4,5,7);
    imagesc(cb);
    axis off; colormap gray; axis image
    title(sprintf('cb=(c-b)>%d',dThr))
    % 
    subfig(4,5,8);
    imagesc(mcb);
    axis off; colormap gray; axis image
    title('mcb')
    % 
    subfig(4,5,9);
    imagesc(bcb);
    axis off; colormap gray; axis image
    title('bcb')
    %    
    subfig(4,5,10);
    imagesc(pcb);
    axis off; colormap gray; axis image
    title('pcb=bcb-mcb')
    %% show detected peak
    subfig(4,5,11);
    lN = length(l);
    if lN>0
        imshow(label2rgb(iml,@jet,'k'),'notruesize');        
    end
    title([num2str(lN) ' peak region(s)'])
    
    %% banner
    if 0    
        subfig(16,1,1);
        set(gcf,'numbertitle','off')
        set(gcf,'name',' N-E-O-V-I-S-I-O-N s.r.o');
        set(gcf,'menubar','none');
        axis off
        t = text(-0.13,0.5,'xwtrace.m - Trace tracking for MSVB 2005');
        set(t,'fontweight','bold');
    end
    drawnow;
end
return