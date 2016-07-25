%NMSUP  Sub-Non-maxima suppression on arbitrary neighbourhood, subpixel
%       detection - ! - Destroys the input image - ! -
%       L = snmsup(I,S) finds local maxima of image I on neighbourhood of size S.
%       It works also for detection in rows or columns if one entry of S is zero.
%       Parameters:
%         I ... Size (M,N). The input image.
%         S ... Size (1,2). Size of the neighbourhood.
%         L ... Size (5,?). The list of the maxima found in subpixel precision and 2nd derivatives.
%               L(:,k) = [i,j,f,dii,djj], where
%                 i, j ... indices of the maximum, found by fitting univariate parabolas to 4-neighbourhood of (i,j)
%                 f ... value I(i,j), by bilinear interpolation
%                 dii, djj ... second derivatives at the point (i,j),
%                                dii = d^2I(i,j)/di^2,
%                                djj = d^2I(i,j)/dj^2.
%                              (For detection in rows/columns, corresponding derivative is NaN).
%
%       MEX implementation. Works only for DOUBLE REAL matrices (!).
%
%       Type this file for demo.
%       See also NMSUP

% Tom Werner, Oct 2001, T.Pajdla Oct 2005


% demo:

% prepare image
I = double(randn(100)>.95);
h = [.5 1 .5]; I = conv2(h,h,I);
I = I.^2;

% subpixel peak detection
E = snmsup(I,[2 2]);
plot(E(2,:),E(1,:),'.','markersize',.1)
