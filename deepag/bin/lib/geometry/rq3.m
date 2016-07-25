%RQ3      Pajdla: Returns a 3x3 upper triangular R and a unitary Q so that X = R*Q
%
%       function [R,Q] = rq3(X)
%
%	X 	= input matrix,
%	Q	= unitary matrix
%	R	= upper triangular matrix
%
%       See also QR.

%       Author:         Tomas Pajdla, Tomas.Pajdla@esat.kuleuven.ac.be 
%                                           pajdla@vision.felk.cvut.cz
%                       05/28/94 ESAT-MI2, KU Leuven
%       Documentation:                            
%       Language:       Matlab 4.1, (c) MathWorks 
%       Last change  : $Id: rq3.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : Ready ($Source: /home/cvs/Matlab/utils/rq3.m,v $)                       
%
function [R,Q] = rq3(X)
 
[iR,iQ] = qr(inv(X));
R = inv(iR);
Q = inv(iQ);
