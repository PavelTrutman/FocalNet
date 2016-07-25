% [L,LL] = u2L(u,C) - perspective back-projection
% 
% u ... 2(3) x n image points, 2 -> 1's augmented
% C ... camera C = P \in R^{3 x 4} or 
%              C.K  ... internal calibration matrix
%              C.E  ... camera euclidean coordinate system transformation
%              C.rp ... [a] parameter of the radial distortion division
%                       model r  = r' * 1/(1+a*r'^2)
%              C.irp ... [a] parametr of the inverse radial distortion
%                             r' = r * (1+a*r^2)
%              C.rp has precedence over C.irp if both are present
% L  ... 4 x 4 x n ray Plucker matrix
% LL ... 4 x 4 x n ray dual Plucker matrix

% (c) T.Pajdla, www.neovision.cz, Oct 2, 2004
function [L,LL] = u2L(u,C)

if size(u,1)<3
    u = [u;ones(1,size(u,2))];
end
if ~any(size(C)-[1 1])
    iK        = inv(C.K);
    x         = iK*u;       
    if exist('C.rp')
        r         = sqrt(sum(x(1:2,:).^2));
        x(1:2,:)  = x(1:2,:)./[[1;1]*(1+C.rp*r.^2)];                      % radial distortion    
    elseif exist('C.irp')
        %
        % solve for r 
        %
        % 0 = a*r^3+r-rr = det([ r   0  1
        %                       -rr  r  1
        %                        0 -1/a r])
        % to correct rr as 
        % 
        % eig([ 0   0  1
        %      -rr  0  1
        %       0 -1/a 0])    
        %
        % or by an itterative fixed-point method
        %
        % r_{i+1} = rr_i - a r_i^3
        %
        x0 = x;
        for i=1:1% fixed-point method
            r2        = sum(x(1:2,:).^2);    
            x(1:2,:)  = x0(1:2,:)- x(1:2,:).*[[1;1]*(C.irp*r2)];          
        end
    end
    x         = [x;zeros(1,size(x,2))];   
    x         = inv(C.E)*x;
    C         = inv(C.E)*[0;0;0;1];
    % for i=1:size(x,2)
    %    oL(:,:,i) = XX2L(x(:,i),C); 
    % end
    % A faster implementation of the above
    L = XX2L(x,C); 
else
    for i=1:size(u,2) 
        z        = xx(u(:,i));
        [U,D,V]  = svd(z*C); 
        L(:,:,i) = V(:,3)*V(:,4)'-V(:,4)*V(:,3)';
        LL(:,:,i)= V(:,1)*V(:,2)'-V(:,2)*V(:,1)';        
    end
end

