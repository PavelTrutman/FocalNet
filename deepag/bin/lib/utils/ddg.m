%DDG      Pajdla: Convolves 1-D signal with the 2nd derivative of a gaussian 
%
%	function [ddsignal,filter]=ddg(signal,s)
%
%	signal	= input signal
%	s	= sigma
%	ddsignal= output signal
%	filter	= sampled filter
%	
%       See also:  CONV, CONV2.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			04/12/95 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: ddg.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/ddg.m,v $)  			 
%
function [dim2,filc2]=ddg(im2,s)

 t               = -2*s:2*s;
 fil2            = -2^(1/2)/pi^(1/2)/s^3*exp(-t.^2/s^2)+2*2^(1/2)/pi^(1/2)/s^5.*t.^2.*exp(-t.^2/s^2);
 kor2            = abs( sum((fil2.*(fil2>0)))/sum((fil2.*(fil2<0))));
 filc2           = fil2;
 filc2([fil2<0]) = filc2([fil2<0])*kor2;
 kor             = sum(filc2.*(filc2>0));
 filc2           = 0.5*filc2/kor;
 dim2            = conv2(im2,-filc2,'same');

return
