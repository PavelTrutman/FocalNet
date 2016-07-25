% h = imagescrgb(im) - Imagesc for RGB image 
% 
% im ... image
% h  ... handle to image object

% (c) T.Pajdla
% Aug 9, 2004
function h = imagescrgb(im)

im = double(im);
im = im - min(im(:));
im = im / max(im(:));
h  = imagesc(im);

