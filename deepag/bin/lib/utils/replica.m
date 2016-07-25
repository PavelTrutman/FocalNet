%REPLICA  Pajdla: Replicates a matrix according to the pattern 
%
%
% 	function rv = replica(v,p)  
%
%	 v	=  input matrix
%	 p	=  replication pattern matrix
%	rv	=  replicated matrix
%	
%       See also:  KRON, CHRFILL.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			03/06/95 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: replica.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/replica.m,v $)  			 
%
function rv = replica(v,p)
 
 rv = kron(p,v);

return
