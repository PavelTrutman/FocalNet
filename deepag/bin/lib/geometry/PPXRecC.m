% c = PPXRecC(P1,X1,P2,X2) - Consistency of two Reconstructions: X1, P1{1}, P1{2} vs X2, P2{1}, P2{2}, the smaller, the better 
%
% P1, P2    = {1:2} of 3 x 4 camera projection matrix
% X1, x2    = 3 x n 3D points
% c         = [c ma1 ma2]
%              c         = consistency abs(ma1-ma2)/min(ma1,ma2)
%              ma1 , ma2 = mean apical angles

% T. Pajdla, pajdla@cmp.felk.cvut.cz, 2015-09-05
function c = PPXRecC(P1,X1,P2,X2)

a1 = PPX2ae(P1{1},P1{2},X1); % apical angles
if isempty(a1), ma1 = NaN; else ma1 = mean(a1(3,:)); end % means at points only
a2 = PPX2ae(P2{1},P2{2},X2); % apical angles
if isempty(a2), ma2 = NaN; else ma2 = mean(a2(3,:)); end % means at points only
if all(isfinite([ma1 ma2]))
    c  = [abs(ma1-ma2)/min(ma1,ma2) ma1 ma2];
else
    c  = [Inf ma1 ma2];
end