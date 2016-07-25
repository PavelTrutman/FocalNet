%CORTON   Pajdla: Closest ortonormal matrix with one defined row
%
%	function [S,fe] = corton(R,n)
%       
%       R   = 3-D rotation matrix
%       n   = 3-D vector
%       S   = rotation matrix, S = [s1 s2 n]
%             norm(R(:,1:2)-S(:,1:2),'frobenius') -> min for all S \in SO(3)
%       fe  = frobenius error fe = norm(R(:,1:2)-S(:,1:2),'frobenius')
%
%
%       See also:  SVD, doc/cortog3.ma

%	Author: 	Tomas Pajdla, pajdla@vision.felk.cvut.cz
%			04/5/95 Computer Vision Laboratory, CTU Prague
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: corton.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/corton.m,v $)  			 
%
function [S,fe]=corton(R,n)

 siz = size(R);
 if(siz(:)~=[3 ; 3])
  error('corton: Error: R must be 3 x 3'); 
 end

 n = n(:);

 if(size(n)~=[3])
  error('corton: Error: n must be 3 x 1'); 
 end

 if(([(R*R')-eye(3)]>2*eps) | ([(R'*R)-eye(3)]>2*eps))
  error('corton: Error: R must be an orthogonal matrix');
 end

 if(norm(n)<2*eps)
  error('corton: Error: n is close to zero'); 
 end
 n = n/norm(n);

 m   = R'*n;
 den = (m(1)^2+m(2)^2);
 U   = [ [(m(2)^2+m(1)^2*m(3))/den,(m(1)*m(2)*(-1+m(3)))/den,-m(1)]',...
         [(m(1)*m(2)*(-1+m(3)))/den,(m(1)^2+m(2)^2*m(3))/den,-m(2)]',...
         m 
       ];
 S   = R*U;
 fe  = norm(R(:,1:2)-S(:,1:2),'fro');
return
