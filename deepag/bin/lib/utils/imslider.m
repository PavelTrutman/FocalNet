% imslider(im,lbl,n) - image slider
%
% im  = {} of images
% lbl = {} of text labels to be ploted under images, can be omitted or []
% n   = the number of visible images, default = 3 

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 3 Nov 2009
function imslider(im,lbl,n)

if nargin<3
    n=3;
end
if nargin<2 || isempty(lbl)
    lbl = num2str((1:length(im))');
    lbl = num2cell(lbl,2);
end
if ~isempty(im)
    set(gcf,'toolbar','figure');
    pos=get(gcf,'position');
    us = uicontrol('Style','slider','Position',[0 0 pos(3) 15]);
    iminf.imN=length(im);
    iminf.imh=size(im{1},1);
    iminf.imw=size(im{1},2);
    iminf.lpos = [0 cumsum(cellfun('size',im,2))];
    iminf.uch=us;
    iminf.n = n;
    set(gcf,'userdata',iminf);
    set(us,'Min',1,'Max',iminf.imN,'Value',1,'Callback',@imsldruchclbk,'SliderStep',[1 n]/iminf.imN);
    imagesc(cat(2,im{:})); axis image; axis off;
    imsldruchclbk();
    for k=1:length(im)
        text(iminf.lpos(k),size(im{k},1)+30,sprintf('%s',lbl{k}));
    end
end
return


