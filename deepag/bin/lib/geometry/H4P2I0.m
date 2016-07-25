% H = H4P2I0(P) - Matrix H transforming 3x4 camer matrix P of a finite cameras  into [I|0]=P*H 
%
% P = 3x4 rank 3 projection matrix
% H = 4x4 matrix H(4,:)=[0 0 0 1] sich that [I|0]=P*H

% T.Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-08-23
function H = H4P2I0(P)

if rank(P(:,1:3))>2
    H = [P(:,1:3)\[eye(3) -P(1:3,4)];[0,0,0,1]];
else
    H = [];
end
