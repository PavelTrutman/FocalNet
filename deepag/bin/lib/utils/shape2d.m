%SHAPE2D Pajdla: Computes a shape from unorganized points in 2D
%                          	                                    		
%	function  e = shape2d(x)
%
% 	x   =  coordinates of vertices of triangulation [x1,y1 ; x2,y2 ; ...]
%	e   =  edges ~ pairs of indices into x	
% 
%	See also: DEALUNAY

%	Author: 	Tomas Pajdla, pajdla@vision.felk.cvut
%
%			01/28/97 Center for Machine Perception, CTU Prague
%	Documentation:                 	 	  
%	Language: 	Matlab 4.1, (c) MathWorks
%       Last change  : $Id: shape2d.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/shape2d.m,v $)  			 
%
function [e,ev,ei,v] = shape2d(x)
 
 v = delaunay(x);
 %
 % find the shortest and the second shortest edge to
 % all the vertices in delaunay triangulation

 e       = [v(:,[1 2]); v(:,[2 3]) ; v(:,[3 1])];
 ev      = [e ; e(:,[2 1])];
 [tmp,i] = sort(ev(:,1)*10^floor(log10(max(ev(:,2)))+1)+ev(:,2));
 e       = ev(i,:);
 i       = ~[0 ;diff(e(:,1))==0 & diff(e(:,2))==0];
 e       = e(i,:);
 xe      = x(e(:,1),:)-x(e(:,2),:);
 d       = sqrt(sum([xe.*xe]')');    
 s       = 10^floor(log10(max(d))+1)*e(:,1)+d;
 [s,i]   = sort(s);
 ei      = e(i,:); 
 i       = [1; diff(ei(:,1))];
 i       = i | [1 ; i(1:length(i)-1)];
 e       = ei(i,:);
return







