%polylinx - Pajdla: Coordinates of points sampled from a polyline
%
%     function x = polylinx(v,d,m)
%
%     v     = [ [x y ... z    ,   w       ]' ; ...] , size(v)=[k+1 n], 
%               [k-dim vertex ,  validity ]' ; ...]
%             w = 1 ... valid
%                 0 ... invalid
%
%               Only valid vertices are interpolated
%
%               v = [ x_1 x_x2 x_3 x_4 x_5 x_6 
%                      1   1    0   1   1   1  ]
%
%                      *---*        *---*---* 
%
%     d     = sampling distance
%     m     = interpolation method
%             'lin' ... linear sampling, in constant distance;
%                       vertices are preserved, last distance 
%                       on each line <= d
%                       (implicit)
%             'bre' ... 4-neighbourhood digital line, Only for a 2D curve !
%     x     = coordinates of the points on the polyline 
%             x = [[x_1; x_2;...;x_k; n_1], ...]
%     n		= normal at each point to the curve projected to the plane of 
%             the first two coordinates
%             n = [[n_1; n_2], ...]
%
% 
%     See also:  . 

%       Author       : Tomas Pajdla, pajdla@cmp.felk.cvut.cz
%                      Apr 4, 1999, Center for Machine Perception, 
%                      Czech Technical University, Prague
%       Documentation:                            
%       Language     : Matlab 5.1, (c) MathWorks                         
%       Last change  : $Id: polylinx.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/polylinx.m,v $)
%
function [x,n] = polylinx(v,d,m,dbflag)

if nargin<4
	dbflag = 0;
end

if nargin<3
	m = 'lin';
end

if nargin==0
   fi= 0:0.1:16*pi;
   v = [fi.*cos(fi);fi.*sin(fi);ones(1,size(fi,2))];
   d = 1;
   m = 'bre';
   dbflag = 0;
end

if size(v,2)<1
   return   
end

w     = v(end,:);
v     = v(1:end-1,:);

if ~any(w)
   return;
end

if strcmp(m,'lin')
   tv = diff(v')';
   sv = sqrt(sum(tv.^2)); 
   tv = tv./(ones(size(tv,1),1)*sv);
   x  = [];
	for i=1:size(v,2)-1
      if w(i)
	      s = 0:d:sv(i);
		   s = s(s<sv(i));
         x = [x [      v(:,i)      *ones(1,size(s,2)) + tv(:,i)*s 
                 [tv(2,i);-tv(1,i)]*ones(1,size(s,2))            ]];
      elseif ~w(i+1)
         x = [x [ v(:,i+1) ; [tv(2,i);-tv(1,i)]]];
      end 
   end
   if w(i+1)
         x = [x [ v(:,i+1) ; [tv(2,i);-tv(1,i)]]];
   end
end

if strcmp(m,'bre')
   if size(v,1)~=2
      error('Error: v is not of dimension 2 in polylinx(...,''bre'')');
   end
   vo  = v;
   % round vertices to a digital grid
   v   = round(v);
   % remove invalid repetitive vertices 
   tv  = diff(v')';
   sv  = sqrt(sum(tv.^2)); 
   ix  = [1 sv>0] | w; 
   v   = v(:,ix);
   w   = w(ix);
   % remove all repetitive vertices => pixel is valid if there was at 
   % least one valid vrtex
   tv  = [diff(v')' NaN*ones(size(v,1),1)];
   sv  = [sqrt(sum(tv.^2))]; 
   ix  = logical([1 sv(1:end-1)>0]); 
   v   = v(:,ix);
   w   = logical(w(ix));
   tv  = [diff(v')' NaN*ones(size(v,1),1)];
   sv  = [sqrt(sum(tv.^2))]; 
   tnv = tv./(ones(size(tv,1),1)*sv);
   % indexes of neighbouring vertices; itself if the neighbour is invalid
   px  = [1 [[2:size(v,2)  ] - 1 + ~w(1:end-1)]          ];
   ax  = [  [[1:size(v,2)-1] + 1 - ~w(2:end  )] size(v,2)];
   % normalized vector to the previous vertex; normalized vector to the nest vertex
   tpx = v-v(:,px);
   tpx = sdiv(tpx,(ones(size(tpx,1),1)*[sqrt(sum(tpx.^2))]));
   tax = v(:,ax)-v;
   tax = sdiv(tax,(ones(size(tax,1),1)*[sqrt(sum(tax.^2))]));
   % normals and mean normal angles at each vertices 
   nv   = [0 1;-1 0]* tnv;
   fnv  = atan2(nv(2,:),nv(1,:));
   fd   = sum([any(tpx);any(tax)]);
   fnvx = (atan2(tpx(2,:),tpx(1,:))+atan2(tax(2,:),tax(1,:)))./fd-pi/2;
   % mark stats of valid segments
   bx   = [1 abs(diff(w))];
   % remove invalid vertices
   vx   = v(:,w);	    % valid vertices
   tx   = tv(:,w);    % tangent vectors
   tnx  = tnv(:,w);   % normalized tangent vectors
   sx   = sv(:,w);    % lengths of the tangent vectors
   fnx  = fnv(:,w);   % normal vectors
   fnvx = fnvx(:,w);  % interpolated normal vectors
   bx   = logical(bx(:,w));    % start of connected segments
   ix   = 1:size(bx,2); 
   lx   = zeros(size(bx));
   lx(bx) = diff([ix(bx) ix(end)+1])';    % lengths of the segments
   if dbflag % debug plots 
      %      figure
      hold on
      nvx  = [cos(fnvx);sin(fnvx)];   
      w    = logical(w(w));
      plot(v(2,:),v(1,:));
      hold on
      plot(vo(2,:),vo(1,:),'r--');
	   plot(vx(2,w),vx(1,w),'.','markersize',20);
	   plot(vx(2,~w),vx(1,~w),'o');
	   quiver(v(2,:),v(1,:),nv(2,:),nv(1,:),0.2,'b');
      quiver(vx(2,:),vx(1,:),nvx(2,:),nvx(1,:),0.3,'r');
      axis equal
   end
   % generate lines
   if 1
   x = [];
   n = [];
   l = [];
   k = 1;
   for i=1:sum(bx)
      for j=1:lx(k)-1
         x = [x vx(:,k)];
         n = [n fnvx(:,k)];
         % Bressenham
         switch sx(k)
         case 1,    % 4-neighbourhood
         case sqrt(2),
            x = [x vx(:,k)+sqrt(2)/2*[1 1;-1 1]*tnx(:,k)];
            n = [n fnx(:,k)];
         case 2,    % 8-neighbourhood
            x = [x vx(:,k)+tnx(:,k)];
            n = [n fnx(:,k)];
         otherwise, % diagonal
            if abs(tnx(1,k))==abs(tnx(2,k))              
               l1 = abs(tx(1,k));
               x1 = repmat([sqrt(2)/2*[1 1;-1 1]*tnx(:,k),...
                     sqrt(2)/2*[1 -1;1 1]*tnx(:,k)],1,l1);
               x1 = cumsum(x1')';
               x  = [x [vx(1,k)+x1(1,1:end-1) ; vx(2,k)+x1(2,1:end-1)]];
               n  = [n fnx(:,k)*ones(1,2*l1-1)];
            else   % other
               if norm(tx(:,k))==norm([1 2])
                  tx(:,k);
               end
               ftx   = rem(atan2(tx(2,k),tx(1,k))+2*pi,2*pi);
               ftx   = fix(ftx/(pi/2))*pi/2;
               Rtx   = [cos(ftx) sin(-ftx); sin(ftx) cos(ftx)];
               rtx   = round(Rtx'*tx(:,k));
               if   rtx(1)>rtx(2)
                  t1 = rtx(1);
                  t2 = rtx(2)/rtx(1);
               else % rtx(1)<rtx(2)
                  t1 = rtx(2);
                  t2 = rtx(1)/rtx(2);
               end
               s     = 0:t1;
               u     = round(t2*s);
               d     = [0 diff(u)];
               si    = s+cumsum(d)+1;
               t     = zeros(1,size(s,2)+sum(d));
               t(si) = d;
               s     = cumsum(~t(2:end-1));
               t     = cumsum( t(2:end-1)); 
               if   rtx(1)>rtx(2)
                  x     = [x [vx(:,k)*ones(1,size(s,2))+Rtx*[s;t]]];
               else % rtx(1)<rtx(2)
                  x     = [x [vx(:,k)*ones(1,size(s,2))+Rtx*[t;s]]];
               end
               n     = [n fnx(:,k)*ones(1,size(s,2))];
            end
         end
         %
         k = k + 1;
      end
      x = [x vx(:,k)];
      n = [n fnvx(:,k)];
      k = k + 1;
   end
   % normalized normal vector
   n  = [cos(n) ; sin(n)]; 
   sn = sqrt(sum(n.^2));
   n  = n./(ones(size(n,1),1)*sn);
   if  dbflag % debug plots 
      plot(x(2,:),x(1,:),'g');
	   plot(x(2,:),x(1,:),'g.','markersize',15);
      quiver(x(2,:),x(1,:),n(2,:),n(1,:),0.3,'g');
      axis equal
      zoom
      grid
   end
end
end

return
