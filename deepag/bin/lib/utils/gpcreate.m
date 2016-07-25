%GPCREATE Pajdla: Creates a Gaussian pyramid
%
% 	p = gpcreate(im,levels,flt)
%
%	Each level of the pyramid p is created by the subsampling of the previous 
%	level convolved with the filter flt.
%
%	im 	= input image
%	flt	= filter matrix, implicitly flt = [1/4 1/2 1/4]' * [1/4 1/2 1/4];
%	levels	= number of levels in the pyramid
%
%	for method see MASKFILT
%
% See also: PCREATE, SPAR, PGETLVL

%	Author       : R. Sara, V. Smutny, T.Pajdla, pajdla@vision.felk.cvut.cz
%                      08/11/96 Computer Vision Laboratory, 
%                      Czech Technical University, Prague
%	Documentation:                 	 	  
%	Language     : Matlab 4.2, (c) MathWorks  			 
%       Last change  : $Id: gpcreate.m,v 1.2 2007/02/26 10:40:56 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/gpcreate.m,v $)
%
function pyr = gpcreate(im,levels,flt)

step = [2 2];

if nargin<3
 flt =  [1/4 1/2 1/4]' * [1/4 1/2 1/4];
end

if step.^levels > size(im)
 error('Too many levels or too large step')
end

if prod(size(step))==1
 step = [step step];
end

imsiz = size(im);
 % pyramid total size:
pyrsiz = 5;	
for l = 1:levels
 pyrsiz = pyrsiz + prod(imsiz);
 imsiz = floor(imsiz./step); 
end

pyr = zeros(pyrsiz,1);
pyr(1:2) = size(im);	% image size
pyr(3) = levels;	% how many levels in the pyramid
pyr(4:5) = step;	% the shrinkage factors between two levels

beginn = 6;
endn = beginn+prod(size(im))-1;
pyr(beginn:endn) = im(:);
for l = 1:levels
 im = conv2(im,flt,'same');
 im = im(1:2:end,1:2:end);
 beginn = endn+1;
 endn = beginn+prod(size(im))-1;
 pyr(beginn:endn) = im(:);
end

return
