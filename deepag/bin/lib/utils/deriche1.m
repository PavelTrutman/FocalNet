function O = deriche1(I,alpha)

% DERICHE1  1-D (row-wise) Deriche edge detection.
%           Y = deriche1(X,alpha) performs convolution of rows of I with 
%           the Deriche filter of impulse response 
%             h(n) = -alpha^3/4*n*exp(-alpha*abs(n)).
%           (Constant -alpha^3/4 is to normalize Integrate[Integrate[h(x)],{1,-Inf,Inf}]] to 1.)
%
%           See also ROWMAX.

[M,N] = size(I);
OP = zeros(M,N);
OM = zeros(M,N);

a = exp(-alpha);
b1 = -2*a;
b2 = a*a;

for n = 3 : N,
  OP(:,n) = I(:,n-1) - b1*OP(:,n-1) - b2*OP(:,n-2);
end
for n = N-2 : -1 : 1,
  OP(:,n) = I(:,n+1) - b1*OP(:,n+1) - b2*OP(:,n+2);
end

for n = N-2 : -1 : 1,
  OM(:,n) = -I(:,n+1) - b1*OM(:,n+1) - b2*OM(:,n+2);
end
for n = 3 : N,
  OM(:,n) = -I(:,n-1) - b1*OM(:,n-1) - b2*OM(:,n-2);
end

O = -alpha^3/4*a*(OP+OM);
