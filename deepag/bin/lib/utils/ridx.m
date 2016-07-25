%RIDX     Pajdla: Non-integer reindexing
%
%	function  oidx = ridx(I,i)
%	
%	I	= m x 2 matrix ~ [ii oi ; ...]
%	i	= input index
%	oidx    = oi(closest(i,ii))
%
%       See also:  VIDX,MATIMAT.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			02/06/95 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: ridx.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/ridx.m,v $)  			 
%
function  oidx = ridx(I,i) 
 [u,ui] = unival(I(:,1));
 I      = I(ui,:);
 i      = i((i>min(I(:,1)))&(i<max(I(:,1))));
 ini    = interp1(I(:,1),I(:,2),i,'linear');
 oidx   = round(ini);
return
