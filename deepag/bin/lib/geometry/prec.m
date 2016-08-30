%PREC    Pajdla: Projective reconstruction from correspondences
%
%       function [X,M1,M2,dx] = prec(u,[K1,K2,E1,E2]) 
%       
%       u	= correspondences in the format [u1 ; v1 ; u2 ; v2]
%	K1,K2	= camera calibration matrices, 3x3 matrices
%	E1,E2	= camera positions, 4x4 matrices repr. an Euclidean
%	          motion 
%	X	= projective reconstruction [x ; y ; z ; w]
%	dx	= the shortest distance between the rays in the Euclidean
%		  reconstruction 
%	M1,M2	= projection matrices, ui = Mi * X
%
%	If [K1,K2,E1,E2] is computed. Otherwise, some projective
%	reconstruction is given. 
%
%       SEE also:  UU2Q, RECONF

%	Author       : Tomas Pajdla, pajdla@vision.felk.cvut.cz
%                      02/25/96 Computer Vision Laboratory, 
%                      Czech Technical University, Prague
%	Documentation:                 	 	  
%	Language     : Matlab 4.2, (c) MathWorks  			 
%       Last change  : $Id: prec.m,v 1.1 2005/04/28 16:54:38 pajdla Exp $
%       Status       : ($Source: /home/cvs/Matlab/geometry/prec.m,v $)
%
function [X,M1,M2,dx] = prec(u,K1,K2,E1,E2)

if nargin>0
  
 if nargin<2
   K1 = eye(3);
 end
 if nargin<3
   K2 = eye(3);
 end
 if nargin<4
   E1 = eye(4);
 end
 if nargin<5
   E2 = eye(4);
 end
 
 %%
  %  Ei = [Ri -Ri*ti
  %        0       1]
  %
  %  X  = [x ; w]
  %
  %  ui = Ki * [Ri -Ri*ti] * X = Ki * Ri * [I -ti] * [x ; w]
  %                            = Ki * Ri * (x - w * ti)
  %
  %  uui = Ri'*inv(Ki)*ui   = x - w * ti
  %
  %%
 
 u1  = e2p(u(1:2,:));	% from affine to homogenous
 u2  = e2p(u(3:4,:));
 R1  = E1(1:3,1:3);
 R2  = E2(1:3,1:3);
 
 if nargin==5		% Euclidean reconstruction
   t1     = E1(1:3,4);
   t2     = E2(1:3,4);
   [x,dx] = reconf(t1,R1,K1,u1(1:2,:),t2,R2,K2,u2(1:2,:))
   M1     = K1*R1*[zeros(3) -t1];
   M2     = K2*R2*[zeros(3) -t2];
   X      = e2p(x);
 else			% projective rec. 
   %%
    % Undo as much calibration as it is known
    %%
   uu1        = R1'*inv(K1)*u1;
   uu2        = R2'*inv(K2)*u1;
   %%
    % compute ti from epipoles   
    %%
   A1         = normu(uu1); 
   nu1        = p2e(A1*e2p(uu1));
   A2         = normu(u2); 
   nu2        = p2e(A2*e2p(uu2));
   nQ         = uu2q('ls',nu1,nu2);
   Q          = A1'*nQ*A2;
   [U1,D1,e1] = svd(Q');
   e1         = e1(:,3);
   [U1,D1,e2] = svd(Q);
   e2         = e2(:,3);
   clear U1 D1 U2 D2 nu1 nu2
   %%
    % Reconstruct
    %
    % M1  = [I   0]
    % M2  = [A -e2]
    %
    % uu1 = [I  0 ] * [x   ; w]  => X ~ [ uu1 ; w ]
    % uu2 ~ [A -e2] * [uu1 ; s]  = A * uu1 - s * e2
    % 
    % [ uu2 ]    [ a1' ]   [ uu1 ]       [ e21 ]
    % [ vv2 ] ~  [ a2' ] * [ vv1 ] - s * [ e22 ]
    % [  1  ]    [ a3' ]   [  1  ]       [ e23 ]
    %
    % (a1'*[uu1 vv1 1]' -s * e21) - uu2 * (a3'*[uu1 vv1 1]' - s * e23) = 0 
    % (a2'*[uu1 vv1 1]' -s * e22) - vv2 * (a3'*[uu1 vv1 1]' - s * e23) = 0
    %
    % a1'*[uu1 vv1 1]'                  - uu2 * a3'*[uu1 vv1 1]' + s * uu2 * (e23 - e21) = 0 
    %                  a2'*[uu1 vv1 1]' - vv2 * a3'*[uu1 vv1 1]' + s * vv2 * (e23 - e22) = 0
    %
    % [uu1 vv1 1           -uu2*uu1 -uu2*vv1 -uu2  uu2*(e23-e21)] [a1]
    % [          uu1 vv2 1 -vv2*uu1 -vv2*vv1 -vv2  vv2*(e23-e22)]*[a2] = 0  
    %                                                             [a3]
    %                                                             [s ]
    %%

    %%
     % Might be better to get A from Q and some other constraints:
     %
     % uu1' * Q * uu2 = 0
     %
     % uu1' * Q * A * uu1 - s * Q * e2 = 0
     %
     % uu1' * Q * A * uu1 = 0 => diag( Q' * uu1 * uu1' * A') = 0
     %
     % diag( G * A' ) = 0 <--- 3 linear equationns for A
     %
     % The rows of A are perp. to the rows of G. The rest of A can be
     % chosen arbitrarily. Can A be a rotation matrix? Q has rank 2
     % iff cameras are translated. In general one wish to have
     % uu1*uu1' of full rank and therefore G is of rank 2. Then a1 can
     % be the right singular vector of G
     % and 
     %
     % 
     
     
    
   uu1 = uu1';
   uu2 = uu2';
   
   C = [uu1              zeros(size(uu1)) -uu2(:,1)*[1 1 1] uu2(:,1)*(e2(3)-e2(1)) 
        zeros(size(uu1)) uu1              -uu2(:,2)*[1 1 1] uu2(:,1)*(e2(3)-e2(2))]; 
   [U,D,a] =           
    
   %% 
    % Compute projection matrices
    %%
   M1 = [eye(3) zeros(3,1)];    
   M2 = [[a1 ; a2 ; a3 ] e2];
   %%
    %  Add known calibration to the projection matrices
    %%
   
   M1 = K1 * R1 * M1;
   M2 = K2 * R2 * M2;
   
 else			

 end
 


else	% demo of the function
  disp('prec: DEMO of the projective reconstruction from 2 views');
end

return




