% A = Hctr(H,x) - Homography mapping covariance trannsformation for point x
%
% H     ... 3 x 3 homography matrix
% x     ... (3)2 x 1 2D point 
% A     ... 2 x 2 transformation matrix
%
% e_x' * inv(C_x) * e_x = I    and    e_y' * inv(C_y) * e_y = I
%
% C_y = T * C_x * T'

% (c) T.Pajdla, www.neovision.cz, Nov 30 2005
function A = Hctr(H,x)
    if size(x,1)==3
        x = h2a(x);
    end
    x0 = H*a2h(x);
    x0 = x0./(ones(3,1)*x0(3,:)).^2;
    A  = [-H(3,1)  0      H(1,1) -H(3,2)  0      H(1,2)
           0      -H(3,1) H(2,1)  0      -H(3,2) H(2,2)];
    A  = A*[x0 zeros(3,1); zeros(3,1) x0];   
return
