% E = XY2rts(X) - E = arg min ||X(4:6,:) - (s*r*X(1:3,:) + t)||
%
% X(1:3,:) ... 3 x n points
% X(4:6,:) ... 3 x n points
% E.r  ... 3 x 3 rotation
% E.t  ... 3 x 1 translation
% E.s  ... 1 x 1 scale

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 6 Jan 2009
function E = XY2rts(X)
X1 = X(1:3,:);
X2 = X(4:6,:);
% centroids
X10 = mean(X1,2);				        
X20 = mean(X2,2);
% centralized points
X1 = X1-X10*ones(1,size(X1,2));
s1 = sum(sqrt(sum(X1.*X1)));
X2 = X2-X20*ones(1,size(X2,2));
s2 = sum(sqrt(sum(X2.*X2)));
% optimal rotation
E.r = XY2rot(X1,X2);
% scale 
E.s = s2/s1;
% translation
E.t = X20 - E.s*E.r*X10;

