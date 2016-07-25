%GAUSMIX Pajdla: Creates a vector of samples from a gaussian PDF mixure
%
%
% function r = gausmix(len,X,S,noiseS)
%
%	len	= length of the row
%	X	= vector of mean values for G. PDF
%	S	= vector of sigma for G. PDF
%	noiseS	= sigma for bacgroung Gaussian noise
%		  additive noise
%
% 	The vector is constructed as a sum of background
% 	gaussian noise with the sigma noiseS and a finite number 
% 	of gaussian distributions each defined by its mean
% 	X(i) and sigma S(i).
%
%       See also:  SMPLGAUS.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			04/12/95 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: gausmix.m,v 1.2 2007/02/26 10:40:56 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/gausmix.m,v $)  			 
%
function r = gausmix(len,X,S,noiseS)

 gnum=length(X);
 if gnum~=length(S)
  disp('Incompatible length of X and S')
  return
 end

 r = zeros(1,len);

 for i=1:gnum,
  r = r + smplgaus(1:len,X(i),S(i));
 end

 r = r + noiseS*randn(1,len);

return
