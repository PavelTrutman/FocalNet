%VIDX     Pajdla: Indexes matrices by a vector 
%	
%       function b = vidx(a,i[,j])
%	
%	a	= input matrix
%	i,j	= index vectors, must have the same size 
%       or
%  i  = index of size m x n, m ... dimension of a
%                            n ... length of the index
%	b       = a(i,j)
%          or
%  b       = a(i(1,:),i(2,:),...,i(m,:))
%
%       See also  RIDX,MATIMAT.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%				            pajdla@vision.felk.cvut.cz
%			12/04/99 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.1, (c) MathWorks
%       Last change  : $Id: vidx.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/vidx.m,v $)  			 
%
function b = vidx(a,i,j)

if nargin==3
   b = a((j-1)*size(a,1)+i);
else
   cp     = cumprod(size(a));
   cp     = [1 cp(1:end-1)];
   j      = i-1;
   j(1,:) = j(1,:)+1;
   j      = cp*j;
   b      = a(j);
end
return
