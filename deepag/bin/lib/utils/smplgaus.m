%SMPLGAUS Pajdla: Values from the Gaussian distribution function
%
% function  g = smplgaus(x,mu,s)
%
%	x	= abscisae vector
%	mu	= mean value
%	s	= sigma
%
%       See also: GAUSMIX, RANDN.

%       Author:         Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%                                           pajdla@vision.felk.cvut.cz
%                       04/08/95 ESAT-MI2, KU Leuven
%       Documentation:                            
%       Language:       Matlab 4.2, (c) MathWorks
%       Last change  : $Id: smplgaus.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/smplgaus.m,v $)                        
%
function g = smplgaus(x,mu,s)
 
g = (1/sqrt(2*pi)/s * exp(-(x-mu).^2/2/s^2));

return
