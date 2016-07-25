%SDIV     Pajdla: Element-wise division with 1/0 replaced by 0
%
%	function d = sdiv(x,y)
%
% 		d = x ./ y  <=> y != 0
%   	  	  = 0           otherwise
% 
%       See also ./ 

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be
%					    pajdla@vision.felk.cvut.cz
%			04/09/94 ESAT-MI2, KU Leuven 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: sdiv.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/sdiv.m,v $)  			 
%
function d = sdiv(x,y) 

 zero = abs(y)<eps;
 d    = (~zero).*x ./ (y + zero);

return
