% imn = imnorm(im) - Normalize im for imagesc 
% 
% im ... image
% imn ... normalized 0<=imn<=1, double

% (c) T.Pajdla
% June 30, 2007
function im = imnorm(im)

im = double(im);
im = im - min(im(:));
im = im / max(im(:));


