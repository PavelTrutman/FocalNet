% [mx,nx,x,A,o,l,il] = trckspot(im,decim,dthrprc,minR,maxR,dThr,doShow) - detect bright spots
%
% im     ... m x n x 1 image
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
function [mx,nx,x,A,o,l] = trckspot(im,decim,dthrprc,minR,maxR,dThr,doShow)

debugLevel   = 0;
subpixMethod = 'winsormean';

% generate the mask
[c,r] = meshgrid([-2*maxR:-maxR maxR:2*maxR],[-2*maxR:2*maxR]);
[c2,r2] = meshgrid([-maxR+1:maxR-1],[-2*maxR:-maxR maxR:2*maxR]);
msk = [r(:)' r2(:)'; c(:)' c2(:)'];
% detect the spot
[o,A,x] = dtctspot(im,decim,dthrprc,minR,maxR,dThr,msk,doShow);
l = cat(2,o{:});
if isempty(o)
    l =[];
    mx=[];
    nx=[];
else
    l = unique(l(4,:));
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
            ix = sub2ind(size(im),ox(1,:),ox(2,:)); % linear index
            vx = log(double(im(ix)))'; %log of a gaussian is a paraboloid
            mx = [sum((vx*[1 1]).*ox(1:2,:)')/sum(vx)]';
        case 'nonsatmean' % mean of nonsaturated pixels
            ox = o{mi}; % pixels above the threshold in the cutout c.s.
            oxo = ox;
            ix = ox(3,:)<0.90*intmax('uint8'); % non-saturated pixels
            ox = ox(:,ix);
            mx = mean(ox(1:2,:)')';            
        case 'nonsatwmean' % mean of nonsaturated pixels
            ox = o{mi}; % pixels above the threshold in the cutout c.s.
            oxo = ox;
            ix = ox(3,:)<0.95*intmax('uint8'); % non-saturated pixels
            ox = ox(:,ix);
            ix = sub2ind(size(im),ox(1,:),ox(2,:)); % linear index
            vx = log(double(im(ix)))'; %log of a gaussian is a paraboloid
            mx = [sum((vx*[1 1]).*ox(1:2,:)',1)/sum(vx)]';                        
        case 'winsormean' % mean of winsorized pixels
            oxo = o{mi};
            rc = round(nx([1 1 2 2])'+ [-maxR maxR -maxR maxR]);
            ic = im(rc(1):rc(2),rc(3):rc(4));
            ic = conv2(ic,ones(7)/49,'same');
            %h  = hist(ic(:),0:255);
            %ch = cumsum(h);
            %ch = ch/max(ch);
            %b  = [find(diff(ch>0.1)) find(diff(ch<0.9))];
            if dthrprc<1
                dthrprc = dthrprc*255;
            end
            mxv = max(ic(:));
            mnv = min(ic(:));
            mxt = 0.7*mxv;
            mnt = max(1.5*mnv,0.3*mxv);
            it  = ic>mnt & ic<mxt;
            [r,c] = find(it);
            r = r + rc(1)-1;
            c = c + rc(3)-1;
            ox = [r';c'];
            mx = mean([r(:) c(:)])';                        
            if 0 % analyse the sensitivity to the threshold
                x=NaN*ones(256,256);
                y=x;
                for i=1:255
                    for j=i+1:255
                        it = ic>i & ic<j;
                        [r,c] = find(it);
                        if ~isempty(r)
                            x(i,j)=mean(r(:));
                            y(i,j)=mean(c(:));
                        end
                    end
                end
            end
    end
    if 1 %norm(mx-nx,2)>0.5
        subfig(2,2,1);
        imshow(log(im)/log(255),'notruesize');
        hold on
        plot(oxo(2,:),oxo(1,:),'.');        
        plot(ox(2,:),ox(1,:),'.c');
        plot(mx(2,:),mx(1,:),'.g');
        plot(nx(2,:),nx(1,:),'.r');
        axis(rc([3 4 1 2]));
        title(sprintf('[%.0f < imt < %.0f]',mnt,mxt));
        pause(0.1);
        close
    end
end

function mx=qmaxfit(u,v,w)
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
            N  = [2*p(1)    p(2)
                    p(2)  2*p(3)];
            V  = -p(4:5);
            mx = N\V
return