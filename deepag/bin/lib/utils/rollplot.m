%ROLLPLOT Pajdla: PLots vectors stored in rows into the sequence of subplots
%
%
%       function  f = rollplot(x,y,subr,plotstr,label,postCmd)
%
%	x	= idependent variable
%	y	= dependent variable
%	subr	= number of rows in the plot
%	plotstr = for plotstr see PLOT.
%	label	= labels of plots stored in rows, if missing plots are numbered.
%	postCmd = string with command performed for every subplot after plotting
%
%       See also PLOT.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%			08/03/94 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.1, (c) MathWorks
%       Last change  : $Id: rollplot.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/rollplot.m,v $)  			 
%
  function  f = rollplot(x,y,subr,plotstr,label,postCmd)

  f = gcf;

  if size(x)~=size(y)
    error('rollplot:Error - incompatible sizes of input matrices');
  end

  [m,n]=size(x);

  if nargin<5
   label=vec2lbls(1:n);
  end

  if nargin<4
   plotstr='-';
  end

  subc   = ceil(n/subr);
  subIdx = reshape((reshape(1:subr*subc,subc,subr))',subr*subc,1);
  for i = 1:n
   subplot(subr,subc,subIdx(i));
   plot(x(:,i),y(:,i),plotstr);
   if ~isempty(label)
     title(label(i,:));
     grid;
   end 
  if exist('postCmd')==1
   eval(postCmd);
  end
  end
