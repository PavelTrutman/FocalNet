% c = r2c(R) - rotation matrix to Cayley representation
%
% R = 3 x 3 rotation matrix (rotation matrix in a 4 x 4 hom transformation) 
% c = Caylay representation [c1,c2,c3]

% Tomas Pajdla, pajdla@cmp.felk.cvut.cz
% 2016-03-22
function c = r2c(R)
R = R(1:3,1:3);
C = (eye(3)-R)*inv(eye(3)+R);
c = [C(3,2);C(1,3);C(2,1)];