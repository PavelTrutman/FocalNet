function [Mr] = coef_mat_two_focal(x,xp)

if nargin == 0
    x = rand(7,2);
    xp = rand(7,2);
end
symbolic = 1;
if symbolic
    syms w1 w2
    syms tx ty tz
    syms b c d
    
    syms x1 x2
    syms xp1 xp2
    
    x = [x1; x2; 1]
    xp = [xp1; xp2; 1]
    
    K1 = [1 0 0; 0 1 0; 0 0 w1]
    K2 = [1 0 0; 0 1 0; 0 0 w2]
    
    R = [1+b^2-c^2-d^2 2*b*c-2*d 2*b*d+2*c; 2*b*c+2*d 1-b^2+c^2-d^2 2*c*d-2*b; 2*b*d-2*c 2*c*d+2*b 1-b^2-c^2+d^2]
    T = [0 -tz ty; tz 0 -tx; -ty tx 0]
    A = transpose(x)*K1*R*T*K2*xp
    
    [c2,t2] = coeffs(A,[b,c,d,tx,ty,tz,w1,w2])
    
    %%
    %computed c and t
    c2 = [ -xp2, x2, xp1, x1, - x1*xp2 - x2*xp1, -2*x1, 2*x2, 2*x1*xp1 - 2*x2*xp2, 2*x1*xp2, 2, -2*x1*xp1, -2*xp2, -2, -2*x2*xp2, 2*x2*xp1, 2*xp1, -xp2, -x2, xp1, -x1, x1*xp2 + x2*xp1, -2, 2*x2*xp2, -2*x2*xp1, 2*xp1, 2*x1*xp2, -2, -2*x1*xp1, 2*xp2, xp2, x2, -xp1, -x1, x1*xp2 - x2*xp1, 2*x1, 2*x2, - 2*x1*xp1 - 2*x2*xp2, xp2, -x2, -xp1, x1, x2*xp1 - x1*xp2]
 
 
    t2 = [ b^2*tx*w1, b^2*tx*w2, b^2*ty*w1, b^2*ty*w2, b^2*tz, b*c*tx*w2, b*c*ty*w2, b*c*tz, b*d*tx, b*d*ty*w1*w2, b*d*ty, b*d*tz*w1, b*tx*w1*w2, b*tx, b*ty, b*tz*w1, c^2*tx*w1, c^2*tx*w2, c^2*ty*w1, c^2*ty*w2, c^2*tz, c*d*tx*w1*w2, c*d*tx, c*d*ty, c*d*tz*w1, c*tx, c*ty*w1*w2, c*ty, c*tz*w1, d^2*tx*w1, d^2*tx*w2, d^2*ty*w1, d^2*ty*w2, d^2*tz, d*tx*w2, d*ty*w2, d*tz, tx*w1, tx*w2, ty*w1, ty*w2, tz]
 
    %%
else
    %monomials
    %t2 = [ b^2*tx*w1, b^2*tx*w2, b^2*ty*w1, b^2*ty*w2, b^2*tz, b*c*tx*w2, b*c*ty*w2, b*c*tz, b*d*tx, b*d*ty*w1*w2, b*d*ty, b*d*tz*w1, b*tx*w1*w2, b*tx, b*ty, b*tz*w1, c^2*tx*w1, c^2*tx*w2, c^2*ty*w1, c^2*ty*w2, c^2*tz, c*d*tx*w1*w2, c*d*tx, c*d*ty, c*d*tz*w1, c*tx, c*ty*w1*w2, c*ty, c*tz*w1, d^2*tx*w1, d^2*tx*w2, d^2*ty*w1, d^2*ty*w2, d^2*tz, d*tx*w2, d*ty*w2, d*tz, tx*w1, tx*w2, ty*w1, ty*w2, tz]
 
    M = [ -xp(:,2), x(:,2), xp(:,1), x(:,1), - x(:,1).*xp(:,2) - x(:,2).*xp(:,1), -2*x(:,1), 2*x(:,2), 2*x(:,1).*xp(:,1) - 2*x(:,2).*xp(:,2), 2*x(:,1).*xp(:,2), 2*ones(size(x,1),1), -2*x(:,1).*xp(:,1), -2*xp(:,2), -2*ones(size(x,1),1), -2*x(:,2).*xp(:,2), 2*x(:,2).*xp(:,1), 2*xp(:,1), -xp(:,2), -x(:,2), xp(:,1), -x(:,1), x(:,1).*xp(:,2) + x(:,2).*xp(:,1), -2*ones(size(x,1),1), 2*x(:,2).*xp(:,2), -2*x(:,2).*xp(:,1), 2*xp(:,1), 2*x(:,1).*xp(:,2), -2*ones(size(x,1),1), -2*x(:,1).*xp(:,1), 2*xp(:,2), xp(:,2), x(:,2), -xp(:,1), -x(:,1), x(:,1).*xp(:,2) - x(:,2).*xp(:,1), 2*x(:,1), 2*x(:,2), - 2*x(:,1).*xp(:,1) - 2*x(:,2).*xp(:,2), xp(:,2), -x(:,2), -xp(:,1), x(:,1), x(:,2).*xp(:,1) - x(:,1).*xp(:,2)];
    Mr = rref(M);
    
    show_Mr = 0;
    if show_Mr
        figure
        spy(Mr);
    end

end
    
    
    
