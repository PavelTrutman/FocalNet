%VEC2LBLS Pajdla: Converts vector elements into a string matrix
%
%
% function  lbl = vec2lbls(vec,format)
%
%	vec    = vector of numbers ( integers )
%	lbl    = matrix having strings of the same length stored in rows.
%	         It can be combined with the imgtext(_,_,lbl) to label 
%	         pixels in an image.
%	format = see SPRINTF, only one '%` is allowed,
%                if omitted '%.0f' is used.
%
%       See also:  TEXT, IMGTEXT.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			03/06/95 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: vec2lbls.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/vec2lbls.m,v $)  			 
%
function  lbl = vec2lbls(vec,format)
  
 if isempty(vec)
  lbl = [];
  return
 end

 if nargin<2
  format = '%.0f';
 end
 fl      = length(format);
 textend = findstr(format,'%');
 if isempty(textend)
   error('vec2lbls:Error - invalid format - character `%` has not been found.');
 end
 textend = textend - 1;

 if format(fl-1)=='%'
  format = ['%.4' format(fl)];
 end
 tmps1 = sprintf(format(textend+1:fl),max(fix(vec)));
 tmps2 = sprintf(format(textend+1:fl),min(fix(vec)));
 tmp  = [length(tmps1);length(tmps2)];
 slen = max(tmp);

 evals=(['lbl =  sprintf(''' format(1:textend) '%' sprintf('%d',slen) format(fl-2:fl) ''',vec(:));']);
 eval(evals);
 lbl  = (reshape(lbl,slen+textend,length(vec)))';

return
