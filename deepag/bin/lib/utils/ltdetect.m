%ltdetect - Laser trace subpixel detection (columnwise). 
%
%       function [sx,ix,nx] = ltdetect(im,par,iv,xs)
%
%	im		= image 
%	par	= [bt dbt dp df fi]
%	
%		bt		= belief threshold, only maximas > bt are valid
%             if bt<0, returns all ignoring all other conditions
%     dbt 	= minimal difference in belief beween the maximum 
%             on the interval [-dp:0:dp] (i.e. at 0) and on the 
%             borders of the interval (i.e. at -dp and dp).
%             Only maximas foir which im(0)-im(dp)>dbt & im(0)-im(-dp)>dbt
%             are valid.
%		dp 	= 2*dp+1 ... the size of interval on which the parabola 
%			   	   	    is fit to the logarithm of image
%		df 	= 2*df+1 ... the size of interval where the Gaussian
%           	          is computed for filtering of columns, 3*sigma=df
%	           0      ... no filtration 
%     fi = angle in which 1D maximua detection is done [in radians], 
%          [0 pi/4 pi/2 3*pi/2] directions are implemented, other directions
%                               are rounded to the implemented ones
%          column direction (i.e. 0 radans) is implicit
%	iv = image mask size(iv) = size(im), zero parts of iv are ignored
%       in im
%  xs = starting points with angle of normal [[row; col; fi],...]
%	ix	= [row; col; brightness of the detected pixel] 
%                 trace detected in pixel accuracy by finding a pixel with 
%                 the maximal brightness in each column
%
%	sx	= [row ; col; brightness maximum on the fitted parabola; 
%                 trace detected in a subpixel accuracy by fitting a parabola  
%                 in the neighbourhood [-dp:dp] for each ix
%
%  nx = [[nr;nc], ...] estimate of the normal at point sx
%
%       See also:  . 

%       Author       : Tomas Pajdla, pajdla@cmp.felk.cvut.cz
%                      Mar 20, 1999 Center for Machine Perception, 
%                      Czech Technical University, Prague
%       Documentation:                            
%       Language     : Matlab 5.1, (c) MathWorks                         
%       Last change  : $Id: ltdetect.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/ltdetect.m,v $)
%
function [sx,ix,nx] = ltdetect(ix,par,iv,xs)

if nargin<2
   bt    = 0;
   dbt   = 0;
   dp    = 3;	
   df    = 0;     
   fi    = 0;
else
   bt    = par(1);
   dbt   = par(2);
   dp    = par(3);
   df    = par(4);
   fi    = par(5);
end

% filter by a gaussian in the column direction to increase smoothness
if df
 sf = df/3; 
 xf = [-df:df];
 f  = (1/sqrt(2*pi)/sf * exp(-(xf.^2)/2/sf^2));
 f  = f/sum(f);
 ix = conv2(ix,f,'same');
end

% get image logarithm to make parabolas form gaussians
ix = log(ix.*(ix>=0)+1);

% mask out unwanted part of the image
if exist('iv')==1
   iv      = iv==0;
   iv(iv)  = -Inf;   
   iv      = ix + iv;
else
   iv      = ix;
end

if 0
   fi      = mod(fi+2*pi,2*pi);
   fi      = rem(round(fi/(pi/4)),4); %[0 pi/4 pi/2 3*pi/4 ...] -> [0 1 2 3]
   % find maximas in direction fi
	if fi==0 % in columns
	   [m,mi] = max(iv);                      % Find maximas
	   mj     = 1:size(iv,2);
	end
	if fi==2 % in rows
	   mi     = 1:size(iv,1);                 % Find maximas
	   [m,mj] = max(iv');
	end
	if fi==1 | fi==3   
	   if fi==1, % in diagonal LU-RD
	      idx0 = [0: 1: size(iv,2)-1 ; 0:size(iv,2)-1];
	      mini = -size(iv,2)+2;
	      maxi =  size(iv,1);
	   end
	   if fi==3  % in diagonal RU-LD
	      idx0 = [0:-1:-size(iv,2)+1 ; 0:size(iv,2)-1];
	      mini = 1;
	      maxi = size(iv,1)+size(iv,2)-1;
	   end
	   j = 1;
	   for i=mini:maxi                                     % generate shift
	      idx = [i ; 1]*ones(1,size(idx0,2)) + idx0;       % shifted diagonal
	      idx = idx(:,idx(1,:)>=1 & idx(1,:)<=size(iv,1)); % clip diagonal
	      vdx = vidx(iv,idx(1,:),idx(2,:));                % sample image on clipped diagonals
	      [mx,mix] = max(vdx);                             % get maximum
	      m(j)     = mx;
	      mi(j)    = idx(1,mix);
	      mj(j)    = idx(2,mix);
	      j        = j+1;
      end
   end
   % Detect the maximum value to a subpixel resolution by fitting a 
	% 1D parabola in a given direction
	
	% generate intervals
	if fi==0
	   x    = [-dp:dp]'*ones(1,size(mi,2)); 
	   y    = zeros(2*dp+1,size(mj,2));
	end
	if fi==2
	   x    = zeros(2*dp+1,size(mi,2)); 
	   y    = [-dp:dp]'*ones(1,size(mj,2));
	end
	if fi==1
	   x    = [-dp:dp]'*ones(1,size(mi,2)); 
	   y    = [-dp:dp]'*ones(1,size(mj,2));
	end
	if fi==3
	   x    = [dp:-1:-dp]'*ones(1,size(mi,2)); 
	   y    = [-dp:dp]'*ones(1,size(mj,2));
	end
	
	xi         = x + ones(2*dp+1,1)*mi;
	xj         = y + ones(2*dp+1,1)*mj;
	
	% clip intervals by image borders
	ixi        = all(xi>=1) & all(xi<=size(iv,1)) & all(xj>=1) & all(xj<=size(iv,2));
	xi         = xi(:,ixi);
	xj         = xj(:,ixi);
	
	% Subpixel detection = fit a parabola in direction fi on the interval [-dp dp]
	dx = x(2,1)-x(1,1);
	dy = y(2,1)-y(1,1);
	sc = sqrt(dx^2+dy^2); % scale of the direction
	
	for i=1:size(xi,2)
	 t(:,i)   = [-dp:dp]'* sc;
	 z(:,i)   = vidx(ix,xi(:,i),xj(:,i));
	 p(:,i)   = polyfit(t(:,i),z(:,i),2)';
	end
	
	ss    = -p(2,:)./p(1,:)/2;
	si    = [ss * dx + xi(dp+1,:) ; ss * dy + xj(dp+1,:)];
	
	% evaluate the belief of the detected positions
	bf  = exp(z(dp+1,:));
	bfm = exp(z([1 end],:));
	
	% remove those which are under the threshold & are not concave & are close to limits
	% & the variatin inside of the interval is high 
	bi = bf>bt & p(1,:)<0 & abs(ss)<=dp-1 & bf-bfm(1,:)>dbt & bf-bfm(2,:)>dbt; 
	ix = [mi ; mj ; exp(m)];
	sx = [si(:,bi) ; bf(bi)];
end

if 1
   % fit 2-D parabolas, and then extract the ridge
   %
   % iv ... masked log of image
   
   if exist('xs')==0
      % find maximas in direction fi
      %
      fi        = mod(fi+2*pi,2*pi);
      fi        = rem(round(fi/(pi/4)),4); %[0 pi/4 pi/2 3*pi/4 ...] -> [0 1 2 3]
		if fi==0 % in columns
		   [m,mi] = max(iv);                      % Find maximas
		   mj     = 1:size(iv,2);
		end
		if fi==2 % in rows
		   mi     = 1:size(iv,1);                 % Find maximas
		   [m,mj] = max(iv');
		end
		if fi==1 | fi==3   
		   if fi==1, % in diagonal LU-RD
		      idx0 = [0: 1: size(iv,2)-1 ; 0:size(iv,2)-1];
		      mini = -size(iv,2)+2;
		      maxi =  size(iv,1);
		   end
		   if fi==3  % in diagonal RU-LD
		      idx0 = [0:-1:-size(iv,2)+1 ; 0:size(iv,2)-1];
		      mini = 1;
		      maxi = size(iv,1)+size(iv,2)-1;
		   end
		   j = 1;
		   for i=mini:maxi                                     % generate shift
		      idx = [i ; 1]*ones(1,size(idx0,2)) + idx0;       % shifted diagonal
		      idx = idx(:,idx(1,:)>=1 & idx(1,:)<=size(iv,1)); % clip diagonal
		      vdx = vidx(iv,idx(1,:),idx(2,:));                % sample image on clipped diagonals
		      [mx,mix] = max(vdx);                             % get maximum
		      m(j)     = mx;
		      mi(j)    = idx(1,mix);
		      mj(j)    = idx(2,mix);
		      j        = j+1;
	      end
      end 
		% mark weak maxima (m - maximum value, [mi mj] - maximum position)
		wi = exp(m)>bt; 
	   % generate dense sample points on a linear approximation of the vertices
	   % and genrate normal directions
		[xs,n] = polylinx([mi ; mj ; wi],1,'bre',0);
   end
   % fit 2D parabolas
   [Pr,Pc] = meshgrid(-dp:dp,-dp:dp);
	pr      = Pr(:);
	pc      = Pc(:);
   PAi     = pinv([pr.^2 pr.*pc pc.^2 pr pc ones(size(pr))]);
   for i=1:size(xs,2)
      % clip the interval
      xi = xs(1,i)+pr; 
      xj = xs(2,i)+pc;
      if all(xi>=1) & all(xi<=size(iv,1)) & all(xj>=1) & all(xj<=size(iv,2))
	      z  = vidx(iv,uint16(xi),uint16(xj));
	      % LSQ fit of a 2D polynomial
	   	P(:,i) = PAi * z;									
	      % 1D polynomial in the direction [u,v] perpendicular to the trace
	   	% z(x,y) = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f
		   % P(:,i) = [a b c d e f]
		   % [u v]  = [pi(3,i)  pi(4,i)]
		   % z(t)   = (v^2*c+a*u^2+b*u*v)*t^2 +(v*e+d*u)*t + f
	      % z(t)   =         p          *t^2 +    q    *t + s
	   	p(:,i)   = [P(3,i)*n(2,i)^2+P(1,i)*n(1,i)^2+P(2,i)*n(1,i)*n(2,i)
		               P(5,i)*n(2,i)+P(4,i)*n(1,i)
	                  P(6,i)];
	      % subpixel correction
	      ss(i)    = -p(2,i)./p(1,i)/2;
	   	si(:,i)  =  ss(i).*n(:,i)+xs(:,i);
	      % the belief of the detected positions
			bf(i)    = exp(polyval(p(:,i),ss(i)));
         bfm(:,i) = exp(polyval(p(:,i),ss(i)+[-dp dp]))';
      end
    end
	% remove those which are under the threshold & are not concave & are close to limits
	% & the variatin inside of the interval is high 
	ix = [mi ; mj ; exp(m)];
   if exist('p')   
      if bt>0
		   bi = bf>bt & p(1,:)<0 & abs(ss)<=dp-1 & bf-bfm(1,:)>dbt & bf-bfm(2,:)>dbt; 
         sx = [si(:,bi) ; bf(bi)];
	 nx = n(:,bi);
      else
         sx = [si ; bf];
	 nx = n;
      end
   end
end
return


