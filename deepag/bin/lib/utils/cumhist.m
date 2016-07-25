%CUMHIST  Pajdla: Plots a histogram and a cummulative histogram
%
%       function [N,X,cX,f]=cumhist(mat,hnum,mode)
%
%	mat	= 1-D or 2-D matrix
%	hnum	= number of bins
%	mode    = 0 ... normalize histogram so that the maximum = 1 (implicit)
%		      1 ... normalize histogram so that the area = 1
%             2 ... do not normalize
%	N,X	= see HIST
%	f	= figure handle
%
%	Histogram is plotted using BAR in yellow scaled to have
%	maximum == 1. Cummulative histogram is scaled to give
%	true Sample Probability Distribution Function.
%
%       See also HIST, HIST2,  HISTO, IMHIST.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%			08/08/94 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.1, (c) MathWorks
%       Last change  : $Id: cumhist.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/cumhist.m,v $)  			 
%
%
function [N,X,cX,f]=cumhist(mat,hnum,mode)
 
 if nargin < 3
     mode = 0;
 end

 if nargin < 2
     hnum = max(20,length(mat)/20);
 end
 
 f           = gcf;
 [N,X] 	     = hist(mat,hnum);
 cX          = cumsum(N);
 hold on
 if mode == 0
     bar(X,N/max(N));
     plot(X,cX/max(N),'-r');
 elseif mode==1
     bar(X,N/cX(end));
     plot(X,cX/cX(end),'-r');
 else
     bar(X,N);
     plot(X,cX,'-r');     
 end
 hold off
  
return
