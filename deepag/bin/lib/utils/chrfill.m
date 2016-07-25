%CHRFILL  Pajdla: Creates and fills the matrix with a given character
%	
%
%       function cm = chrfill(siz,ch)
%
%	siz	= 1 x 2 size matrix
%	ch	= chracter or an integer
%	cm	= string matrix of size "siz" 
%	
%       See also ABS, SETSTR, REPLICA.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be
%					    pajdla@vision.felk.cvut.cz
%			09/02/94 ESAT-MI2, KU Leuven 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: chrfill.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/chrfill.m,v $)  			 
%
function  cm = chrfill(siz,ch) 

 cm = setstr(replica(abs(ch),ones(siz)));

return
