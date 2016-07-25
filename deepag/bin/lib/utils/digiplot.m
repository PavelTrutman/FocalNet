%DIGIPLOT Pajdla: Plots binary digital signal stored in matrix columns 
% 
%
%	function [psig,f]=digiplot(dsig,plotstr)
%
%	dsig	= binary digitaal signal stored on matrix columns
%	plotstr	= plot string, see PLOT
%
%       See also:  PLOT, ROLLPLOT, GEC2BIN, BIN2DEC.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			03/06/95 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: digiplot.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/digiplot.m,v $)  			 
%
function  [psig,f]=digiplot(dsig,plotstr)
 
 [m,n] = size(dsig);
 
 sizes     = (max(dsig)-min(dsig));
 baselines = meshdom(1:2:2*n,1:m);
 baselevel = meshdom(min(dsig),1:m);
 sizes     = meshdom(sizes,1:m);
 
 psig      = (dsig-baselevel)./sizes;
 tmp	   = zeros(size(psig));
 psigf	   = finite(psig);
 tmp(psigf)= psig(psigf);
 psig      = tmp + baselines;
 
 
 if nargout<1
  stairs(psig(:,1));
  hold on
  for i=2:n
    stairs(psig(:,i));
  end
  if nargin==2
   plot(psig,'*g');
  end
  hold off 
  grid;
  axis([1 m 0 2*n]);
 end 

end
