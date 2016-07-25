%PLOTTAB  Pajdla: Plots values from a matrix into a rectangular table
%
%
%       function  f = plottab(a,rLabel,cLabel,hLabel,format,fontSize)
%
%	a	= matrix to be plotted
%	rLabel	= matrix of strings to label rows of the table, 
%		  default = ' '
%	cLabel	= matrix of strings to label columns of the table, i
%		  default = ' '
%	hLabel	= string to label table
%		  default = ' '
%	format	= format for sprintf to plot values in a
%		  default = '%.4f'
%	fontSize= use to set size of of fonts, 
%		  default = 9pt 
%	f	= figure handle
%
%       See also FIGURE, SUBPLOT.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%			08/03/94 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.1, (c) MathWorks
%       Last change  : $Id: plottab.m,v 1.2 2007/02/26 10:40:56 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/plottab.m,v $)  			 
%
function  f = plottab(a,rLabel,cLabel,hLabel,format,fontSize)

if nargin < 2
 rLabel = chrfill([size(a,1) 1],' ');
end

if nargin < 3
 cLabel = chrfill([size(a,2) 1],' ');
end
 
if nargin<4
 hLabel = ' ';
end

if nargin<5
 format = '%.4f';
end

if nargin<6
 fontSize = 9;
end

if ~any(abs(rLabel)-abs(' '))
 cstart = 1;
 rLabel = chrfill([size(a,1) 1],' ');
else
 cstart = 2;
end

if ~any(abs(cLabel)-abs(' '))
 rstart = 1;
 cLabel = chrfill([size(a,2) 1],' ');
else
 rstart = 2;
end

if any(abs(hLabel)-abs(' '))
 cstart = 2; 
 rstart = 2;
end

if size(format,1)==1
 format = chrfill([prod(size(a)) 1],format);
else
 if size(format,1)==size(a,1)
  for i = 1:size(a,1)
   tmp((i-1)*size(a,2)+1:i*size(a,2),:) = chrfill([size(a,2) 1],format(i,:));
  end
  format = tmp;
 else
  if size(format,1)==size(a,2)
   format = chrfill([size(a,1) 1],format);
  end
 end
end

ax      = gca;
delx    = 0.03;
dely    = 0.0;

left    = 0.0;
right   = 1;
bottom  = 0;
top     = 1;
width   = right-left;
height  = top-bottom;

m   = size(a,1)+rstart-1;
n	= size(a,2)+cstart-1;
 
[x,y]   = meshgrid(([1:n]-0.5)*width/(n),([0:m-1]+0.5)*height/(m));
y = flipud(y);

f       = gcf;
axis('off');
 
text(left,y(1,1)+2*dely,hLabel,'units','norm','horiz','left');

if cstart == 2
 for i = rstart:m
  text(x(1,1),y(i,1),rLabel(i-(rstart-1),:),'units','norm','horiz','right','FontSize',fontSize);
 end
end
if rstart == 2
 for i = cstart:n
  text(x(1,i),y(1,1)+2*dely,cLabel(i-(cstart-1),:),'units','norm','horiz','center','FontSize',fontSize);
 end
end

for i = rstart:m
 for j = cstart:n
     if iscell(a) && ischar(a{i-rstart+1,j-cstart+1})
         tx = a{i-rstart+1,j-cstart+1};
     else
         tx   = sprintf(format((i-rstart)*size(a,2)+(j-cstart+1),:),a(i-rstart+1,j-cstart+1));
     end
     text(x(i,j),y(i,j),tx,'units','norm','horiz','center','FontSize',fontSize);
 end
end 

hold('on');
%plot([0 1 1 0 0],[0 0 1 1 0],'-k');

if rstart == 2
% plot([(3*x(1,1)+1*x(1,2))/4 x(1,size(x,2))+x(1,1)],(y(1,1)+y(2,1))/2*[1 1],'-w');
end

if cstart == 2
 if rstart == 1
  upper = (y(1,1)+y(size(y,1),size(y,2))/2);
 else
  upper = (y(1,1)+y(2,1))/2;
 end
% plot((3*x(1,1)+1*x(1,2))/4*[1 1],[upper y(size(y,1),1)-y(size(y,1),size(y,2))/2],'-w');
end
hold('off')

return
