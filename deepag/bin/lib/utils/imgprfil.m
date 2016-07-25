%IMGPRFIL Pajdla: Shows horizontal and vertical profiles in an image
%
%	function imgprfil(m,f1)
%
%	m  	= image matrix
%	f1 	= figure with the image displayed
%
%       See also:  RCPROFIL.

%	Author: 	Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%					    pajdla@vision.felk.cvut.cz
%			04/15/95 ESAT-MI2, KU Leuven
%	Documentation:                 	 	  
%	Language: 	Matlab 4.2, (c) MathWorks
%       Last change  : $Id: imgprfil.m,v 1.2 2007/02/26 10:40:56 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/imgprfil.m,v $)  			 
%
function imgprfil(m,f1)

f1 = subfig(2,2,2,f1);
title('Original image: ');

f2=subfig(2,2,4);
axis('square');
title('Row cut:');

f3=subfig(2,2,1);
axis('square');
title('Column cut:');

n=0;

% Loop, picking up the points.

disp('Input reference point');
but = 1;
while but == 1 
   figure(f1);
   [xi,yi,but] = ginput(1);
   xi = round(xi);
   yi = round(yi);
   n = n + 1;
   if but == 1
     figure(f1);
     title(sprintf('Original image: row=%3.0f, col=%3.0f, val=%.4f',yi,xi,m(yi,xi)));
     
     figure(f2);
     plot(m(yi,:));
     hold on
     plot(xi,m(yi,xi),'+r');
     hold off
     title(sprintf('%s %3.0f','Row cut: row = ',yi));

     figure(f3);
     plot(m(:,xi));
     hold on
     plot(yi,m(yi,xi),'+r');
     hold off
     title(sprintf('%s %3.0f','Column cut: column = ',xi));
   end
end

return
