%POLYAPX  Pajdla: Polygonal approximation of the sampled curve - iterative endpoint fit.
%  
% 	function [seg,err,errs,param] = polyapx(X,merr) 
%
%	X	    = ordered sequence of 3-D points
%	merr	= maximum error allowed
%	vert	= vector of indices to pogon vertices
%	err  	= maximum error 
%	errs	= maximum error for every segment
%	param	= parameterization of data points obtained via their orthogonal
%		      projection onto polygonal segments
%
% See also: 		  

% 	Author: 	Tomas Pajdla, pajdla@vison.felk.cvut.cz 
% 	                08/01/94 ESAT-MI2, KU Leuven 
%
% 	Documentation:  R.O.Duda, P.E.Hart: Pattern Classification and
%			Scene Analysis, John Willey & Sons, 1973.
%			pp. 338-339.
% 
% 	Language: 	Matlab 4.1, (c) MathWorks
%       Last change  : $Id: polyapx.m,v 1.2 2005/09/24 14:37:11 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/polyapx.m,v $)
%                                                            	 
%
  function [vert,err,errs,param] = polyapx(X,merr)
  a          = 0;
  segToCheck = [];
  vert       = [];
  errs		 = [];
  param      = [];

  if ~isempty(X)
   numSeg	      	         = 1;
   numToCheck                = 1;
   segToCheck(numToCheck,:)  = [1 size(X,2)];
   param                     = 0;
  else
   numToCheck      	         = 0;
  end

  while numToCheck > 0
    chS	 = segToCheck(numToCheck,:);
    % the distances of points from the 1st-last-point-line
    V    = X(:,chS(1):chS(2)); % restrict it to the current segment    
    V    = V - V(:,1)*ones(1,chS(2)-chS(1)+1); % move the coordinate sysytem to the first point 
    lv   = V(:,end);           % 1st-last-point-line direction vector 
    lvn  = sqrt((lv'*lv));        % normalize
    if lvn>eps
        lv   = lv/lvn;
    end
    lv   = lv*ones(1,chS(2)-chS(1)+1); % one for each point
    a2   = sum(V.*lv).^2;      % the length of the projection of V onto lv
    c2   = sum(V.*V);          % the length of V
    d    = sqrt(c2-a2);        % the Pythagorean theorem 
    
    [maxDist,maxi] = max(d);
    maxDistIdx	= chS(1) + maxi - 1; 
    a = a+length(d);
    if maxDist > merr+2*sqrt(eps)
      segToCheck(numToCheck,:) 	= [maxDistIdx chS(2)] ;
      numToCheck 		        = numToCheck+1;
      segToCheck(numToCheck,:) 	= [chS(1) maxDistIdx];
    else
      param         =  [param diff(sqrt(a2))];
      vert(numSeg)  =  segToCheck(numToCheck,1);
      errs(numSeg)  =  maxDist;
      numToCheck    =  numToCheck - 1;
      numSeg	    =  numSeg + 1;
    end
   end
  
  vert(numSeg)    =  segToCheck(numToCheck+1,2);
  vert            =  vert';
  errs		      =  errs';
  err             =  max(errs);
  param           =  cumsum(param');





