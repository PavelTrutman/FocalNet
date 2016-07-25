%PFIT     Pajdla: Plane fit using Least Squares or Least Winsorized Squares
%	
%
%       function [p,Vp,rW,res,rs,rres] = pfit(X,W,mode,par)
%	
%	X1		= 3-D points, one point ~ one row
%	W		= a vector of weights, if missingt, ones are used
%	mode		= 0 or missing ... least squares minimization
%					   min sum d^2
%					   par is not used
%			  1            ... least winsorized of squares minimization
%					   par(1) = number of random samples, 
%					            0 means all combinations
%			 	           par(2) = outlier cutoff threshold 
%					   par(3) = robust scale multiplier cut off
%					   r_i is outlier <=>  (r_i > par(3)*robust_scale) || r_i > par(2) 
%					   par(4) = quantile in <0;100>
%							
%			  2 	       ... least winsorized squares followed by 
%					   trimmed total least squares. Outliers are
%					   removed.
%		  
%
%	p		= [nx,ny,nz,d] where [nx ny nz]*[x y z]' + d = 0 
%
%	Vp		= variance matrix of p
%	rs		= robust scale esimate, [1] p. 202, eq. (1.3)
%	rW		= weights after outlier edetection, outliers have weight 0
%	res		= residuals
%	rres    = vector of residuals for the best fit from MLSQ
%	dX		= difference vector
%
%       See also E3, APLE3.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be
%					    pajdla@vision.felk.cvut.cz
%			11/10/94 ESAT-MI2, KU Leuven
%			[1] P.Rousseuw, A.Leroy: Robust Regression and Outlier Detection,
%			  
%         		[2] S. Van Huffel, J. Vanderwalle: The Total Least Squares Problem.
%                	    SIAM Frontiers in Aplied Mathematics 9, ISBN 0-89871-257-0, 1991.		    
%			    John Willey & Sons, 1987, ISBN:0-471-85233-3.
%	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: pfit.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/geometry/pfit.m,v $)  			 
%
function [p,Vp,W,res,rs,rres] = pfit(X,W,mode,par) 
  
 if nargin < 2
  W    = ones(size(X,1),1);
 end
 
 if nargin < 3
  mode = 0;
 end
 
 augOnes = ones(size(X,1),1);
 W       = W./sum(W(:));
 
 %%
  % Less than 3 points are not expected
  %
  %%
 
 if size(X,1) < 3
   res    = 0;
   rs     = 0;
   rres   = 0;
   p      = zeros(4,1);
   Ve     = 'not implemented';
   display('pfit: Warning - less than 3 points supplied');
   return
 end

 %%
  % Robust estimation and outlier detection
  %          Winsorized least squares
  %%

 if (mode == 1)||(mode==2)				  
  samplNum    = par(1);
  cutOff      = par(2);
  rScaleLevel = par(3);
  quantil     = par(4);
  pNum        = size(X,1);			

  rand('seed',sum(100*clock));			% get random point triples of size samplNum, 
						% no two indices in a point triple can be the same
  ridx   = [];
  snKoef = 2;
  while (size(ridx,1)~=samplNum) 
    ridx = round((pNum-1)*(rand(3,snKoef*samplNum)))+1;
    ridx = [ridx ; ridx(1,:)];
    didx = diff(ridx);
    uidx = all(abs(didx));
    ridx = (ridx(1:3,uidx))';
    ridx = ridx(1:min(samplNum,size(ridx,1)),:);
    snKoef = 2 * snKoef;
  end
  
  X1      = X(ridx(:,1),:);				% for all triples
  X2      = X(ridx(:,2),:);
  X3      = X(ridx(:,3),:);

  normal  = cross((X2-X1)',(X3-X1)');
  nnorm   = sqrt(dot(normal,normal));
  nidx    = nnorm>eps;
  normal  = normal(:,nidx);
  nnorm   = nnorm(:,nidx);
  normal  = normal ./ replica(nnorm,ones(3,1));
  d       = - dot(normal,(X1(nidx,:))');
  p       = [normal ; d];
   
  dX      = [X augOnes] * p;
  e       = abs(dX);					% error 
  se      = prctile(e,quantil);

  [mme,mmi] = min(se);					% find minimal winsorized lsq error
  p         = p(:,mmi);
  rres      = dX(:,mmi);
  rs        = 1.4826*(1+5/(max(pNum,7)-6))*sqrt(mme);
  W         = ~((abs(rres)>rScaleLevel*rs)|(abs(rres)>cutOff)).*W;
  W         = sdiv(W,sum(W(:)));
  res       = rres;
 end

								
 %%
  % Least squares solution
  %
  %%

if (mode == 0)||(mode==2)
  if sum(W(:)>eps)>2
   w  = replica(W,ones(1,size(X,2)));
   X0 = sum(w.*X);
   Y  = [X(:,1)-X0(1) X(:,2)-X0(2) X(:,3)-X0(3)];
   
   K       = zeros(3,3);					% direction correlation matrix
   for i=1:size(Y,1)						
     K      = K + W(i)*(Y(i,:))'*Y(i,:);	
   end			

   [v,d]     = eig(K);
   d         = diag(d);
   [md,mdi]  = min(d);
   normal    = v(:,mdi);
   normal    = normal / norm(normal);
   d         = - X0 * normal;
   p         = [normal ; d];
   res       = [X augOnes] * p;
  else
   p   = zeros(4,1);
   res = 0;
  end
 end

 Vp  = 'not implemented';
 
return
