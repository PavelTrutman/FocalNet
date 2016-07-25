%sirpi   - pei: Shift invariant representation of a panoramic image.
%
%       function [ir,fi,fir,e] = sirpi(im [,method,arg1,arg2,arg3])
%
%		im	= image. Must be of odd size. If the size is
%			  even, the last row and column are cutted
%			  off.
%		ir      = shift invariant representation 
%               fi      = final shift of the image
%               fir     = shift of each row
%
%		method	= 'ZPR' ... Image is shifted so that the
%			       	    phase is zero a function computed from rows
%	     
%	                   arg1 = 'mean'   ... of the first harmonics in all rows 
%                                              (equivalent to original ZPR from [1]
%                                 'median' ... of the first harmonics
%                                              in all rows
%			           'wmean' ... weighted mean
%					         from the central row both sides
%					arg3 ... a vector of weights of size(im,1)
%
%		 	           arg2 = vector of row numgers, 
%			                 inf     = all rows ... ZPR on
%				                                 function
%					                         of rows defined by arg1
%	                                 missing = all rows ... original 2D
%		                                                 ZPR uses fft2,
%			                                         equals to:
%				                                 arg1='mean', arg2=inf
%
%				   
%			   e    = energies of of the harmonics of
%				   the function of rows
%
%
%			= 'manual' ... Allows to set the prefered
%			               oriantation into the center of the
%			               by a mouse click

%       See also:  . 

%       Author       : Tomas Pajdla, pajdla@cmp.felk.cvut.cz
%                      Feb 10 1999, Center for Machine Perception, 
%                      Czech Technical University, Prague
%       Documentation:                            
%			@TECHREPORT{Pajdla-TR-170,
%			AUTHOR = {Pajdla, Tom{\'a}{\v s}},
%			INSTITUTION = { FEE CTU },
%			TITLE = {Robot Localization Using Shift Invariant Representation of Panoramic Images},
%			ADDRESS = {FEL {\v C}VUT, Karlovo n{\'a}m{\v e}st{\'\i} 13, Praha, Czech Republic},
%			YEAR = { 1998 },
%			MONTH = { November },
%			NUMBER = { K335-CMP-1998-170 },
%			PSURL = {<a href="ftp://cmp.felk.cvut.cz/pub/cmp/articles/pajdla/Pajdla-TR-170.ps.gz">PostScript</a>},
%			PUBLISHER   = { CMP FEE Czech Technical University },
%			PAGES       = { 23 },
%			}
%
%       Language     : Matlab 5.1, (c) MathWorks                         
%       Last change  : $Id: sirpi.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/sirpi.m,v $)
%
function [ir,fi,fir,e] = sirpi(im,method,arg1,arg2,arg3)

  % initialize arguments
 if nargin<2
   method = 'ZPR';
 end
 if nargin<3
   arg1 = 'mean';
 end
 if nargin<4
   arg2 = NaN; % original ZPR uses fft2
 end  
 if nargin<5
   [m,n] = size(im);   
   arg3  = m-1+mod(m,2);
 end
 
 if strcmp(method,'ZPR')
   if any(isnan(arg2))					
     % Original ZPR
     %
     % make size odd
     [m,n] = size(im);
     im    = im(1:m-1+mod(m,2),1:n-1+mod(n,2)); 
     if any(isinf(arg2))
       arg2 = 1:size(im,1); % ZPR computed from all rows equals to
			    % origianal ZPR but it is computed as a
			    % mean of the first harmonics in all rows
     end  
     % phase of the first harmonics
     f     = fft2(im);
     idf   = -imag(log(f(1,2)));
     % plane of idf slope
     oc    = -floor(size(im,2)/2):floor(size(im,2)/2);
     io    = idf*ones(size(f,1),1)*oc;
     % rotate by multiplying by exp(...*j)  
     g     = ifftshift(exp((io)*j).*fftshift(f));
     % return to spatial domain
     ir    = real(ifft2(g));
     % compute the shift in pixels
     fi    = idf/2/pi*size(im,2);
     % all rows have same shift 
     fir   = fi * ones(size(im,1),1);
     %energies
     e     = NaN * ones(size(im,1),1);
   else
     % ZPR on a function of selected rows
     %
     % make size odd
     [m,n] = size(im);
     im    = im(1:m-1+mod(m,2),1:n-1+mod(n,2));      
     if any(isinf(arg2))
       arg2 = 1:size(im,1); % ZPR computed from all rows equals to
			    % origianal ZPR but it is computed as a
			    % mean of the first harmonics in all rows
     end       
     % select the rows     
     is    = im(arg2,:);
     % get a function of rows
     if strcmp(arg1,'mean') % mean of rows
       mf  = mean(is);
     elseif strcmp(arg1,'median') % median of rows
       mf  = median(is);
     elseif strcmp(arg1,'wmean') % weigted mean of rows
       mf  = mean(is.*(arg3(:)*ones(1,size(is,2))));
     end
     % the first harmonic
     f1    = fft(mf);
     idf   = -imag(log(f1(1,2)));
     f     = fft2(im);
     % plane of idf slope
     oc    = -floor(size(im,2)/2):floor(size(im,2)/2);
     io    = idf*ones(size(f,1),1)*oc;
     % rotate by multiplying by exp(...*j)  
     g     = ifftshift(exp((io)*j).*fftshift(f));
     % return to spatial domain
     ir    = real(ifft2(g));
     % compute the shift in pixels
     fi    = idf/2/pi*size(im,2);
     % all rows have same shift 
     fir   = fi * ones(size(im,1),1);
     % energies
     e     = abs(f1);
   end
 elseif strcmp(method,'manual') % manual alignment
   ir = 64*[im-min(im(:))]/max(im(:)); 
   figure
   plot([size(ir,2)/2 size(ir,2)/2],[1 size(ir,1)],'g');
   hold
   imh = image(ir);
   colormap gray
   axis image
   axis ij
   set(gca,'xticklabel','center')
   set(gca,'yticklabel',[])
   set(gca,'xtick',size(ir,2)/2);   
   axis([1 size(ir,2) 1 size(ir,1)])
   title('left = move to center, right = accept');
   idf = NaN;
   ix  = 1:size(ir,2)-1;
   while 1
     ir  = ir(:,ix);
     delete(imh);
     imh = image(ir);
     drawnow
     ch  = get(gca,'children');
     set(gca,'children',ch([2 1]));
     ylabel(sprintf('shift = %d',idf));
     [c,r,b] = ginput(1);
     if b==3
	 break;
     end
       idf= floor(size(ir,2)/2-c);
       ix = mod([0:size(ir,2)-1]-idf,size(ir,2))+1;
   end
   % compute the shift in pixels
   fi    = idf;
   % all rows have same shift 
   fir   = fi * ones(size(ir,1),1);       
 else 
   error('sirpi: Error - Unknown method!')  
 end     
 return
     
     

