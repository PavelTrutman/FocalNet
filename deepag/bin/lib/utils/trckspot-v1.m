% [mx,nx,x,A,o,l,il] = trckspot(im,decim,dthrprc,minR,maxR,dThr,doShow) - detect bright spots
%
% im     ... m x n x 1 image
% rect   ... valid rectangle in the image [minRow maxRow minCol maxCol]
% decim  ... 1:decim:end image decimation 
% dthrprc... thereshold 0<1  ... percentile 
%                       >=1 ... absolute value 
% minR   ... bright center outer radius
% maxR   ... dark perifery inner radius
% dThr   ... peak decision threshold
% doShow ... 1 = show images
% mx     ... intensity weighted centroid centroid of the object with the maximal area
% nx     ... centroid of the object with the maximal area
% x      ... object centroids
% A      ... object areas
% o      ... object pixel coordinates
% l      ... object labels , 0 = background not included
% il     ... labeled image 

% (c) T.Pajdla, www.neovision.cz, Aug 10, 2005
function [mx,nx,x,A,o,l,il] = trckspot(im,rect,decim,dthrprc,minR,maxR,dThr,doShow)

debugLevel   = 0;
subpixMethod = 'nonsatmean';

[il,l] = dtctspot(im,decim,dthrprc,minR,maxR,dThr,doShow);
for j=1:length(l)
    [r,c] = find(il==l(j));
    ix = sub2ind(size(im),r,c);     
    r  = r + (rect(1)-1);
    c  = c + (rect(3)-1);
    v  = im(ix);
    o{j} = [r';c';v'];
    A(j) = length(r);
    x(:,j) = [mean(r);mean(c)];
end
if isempty(j)
    o=[];
    A=[];
    x=[];
    mx=[];
    nx=[];
else
    [mA,mi] = max(A);
    nx      = x(:,mi);  
    % subpixel by parabola fitting
    if debugLevel>0
        cl = [min(o{mi}(1,:)) max(o{mi}(1,:)) min(o{mi}(2,:)) max(o{mi}(2,:))]; % cutout limits
        ct = NaN*ones(cl(2)-cl(1)+1,cl(4)-cl(3)+1); % undefined cutout
        ox = (o{mi}+1)-cl([1 3])'*ones(1,size(o{mi},2)); % pixels above the threshold in the cutout c.s.
        ix = sub2ind(size(ct),ox(1,:),ox(2,:)); % linear index
        vi = log(im(cl(1):cl(2),cl(3):cl(4))); % log of a gaussian is a paraboloid
        ct(ix) = vi(ix);
        subfig(2,3,1);
        mesh(ct);
    end
    switch subpixMethod
        case 'polynomial' % fitting a plynomial
            ox = o{mi}; % pixels above the threshold in the cutout c.s.
            oxo= ox;
            ix  = sub2ind(size(im),ox(1,:),ox(2,:)); % linear index
            ox  = [ox-nx*ones(1,size(ox,2)); log(double(im(ix)))]'; %zero in maximum, log of a gaussian is a paraboloid
            u  = ox(:,1);
            v  = ox(:,2);
            w  = ox(:,3);
            % Fit a second order polynomial surface
            % w = p(1)*u^2+p(2)*u*v+p(3)*v^2+p(4)*u+p(5)*v+p(6)
            % w = M*p, where
            % p = pinv(M)*w = M\w;
            %
            % Its extreme is obtained by solving:
            % 0 = dw/du = 2*p(1)*u+  p(2)*v+p(4) 
            % 0 = dw/dv =   p(2)*u+2*p(3)*v+p(5) 
            % i.e.
            %  [2*p(1)   p(2)] * [u] = [-p(4)]
            %  [  p(2) 2*p(3)]   [v]   [-p(5)]
            %                N * U   = V
            % a.s.
            % U = N\V  
            % M  = [u.^2 u.*v v.^2 u v ones(size(u))];
            p  = M\w;
            % p  = [p(1);0;p(2:end)];
            N  = [2*p(1)    p(2)
                    p(2)  2*p(3)];
            V  = -p(4:5);
            mx = N\V+nx;
        case 'mean' % the mean value weighted by the intensities
            ox = o{mi}; % pixels above the threshold in the cutout c.s.
            ix = sub2ind(size(im),ox(1,:)-(rect(1)-1),ox(2,:)-(rect(3)-1)); % linear index
            vx = log(double(im(ix)))'; %log of a gaussian is a paraboloid
            mx = [sum((vx*[1 1]).*ox(1:2,:)')/sum(vx)]';
        case 'nonsatmean' % mean of nonsaturated pixels
            ox = o{mi}; % pixels above the threshold in the cutout c.s.
            oxo = ox;
            ix = ox(3,:)<0.90*intmax('uint8'); % non-saturated pixels
            ox = ox(:,ix);
            ix = sub2ind(size(im),ox(1,:)-(rect(1)-1),ox(2,:)-(rect(3)-1)); % linear index
            vx = log(double(im(ix)))'; %log of a gaussian is a paraboloid
            mx = mean(ox(1:2,:)')';            
        case 'nonsatwmean' % mean of nonsaturated pixels
            ox = o{mi}; % pixels above the threshold in the cutout c.s.
            oxo = ox;
            ix = ox(3,:)<0.95*intmax('uint8'); % non-saturated pixels
            ox = ox(:,ix);
            ix = sub2ind(size(im),ox(1,:)-(rect(1)-1),ox(2,:)-(rect(3)-1)); % linear index
            vx = log(double(im(ix)))'; %log of a gaussian is a paraboloid
            mx = [sum((vx*[1 1]).*ox(1:2,:)',1)/sum(vx)]';                        
    end
    if norm(mx-nx,2)>0.5
        subfig(2,2,1);
        imshow(log(im)/log(255),'notruesize');
        hold on
        plot(oxo(2,:),oxo(1,:),'.');        
        plot(ox(2,:),ox(1,:),'.c');
        plot(mx(2,:),mx(1,:),'.g');
        plot(nx(2,:),nx(1,:),'.r');
        axis([nx(2)+[-30 30] nx(1)+[-30 30]]);
        pause(1)
        close
    end
end

