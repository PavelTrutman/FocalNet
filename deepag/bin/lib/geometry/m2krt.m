%[K,R,T] = m2krt(M) - Camera matrix conversion to K,R,T - obsolete chage to P2KRC

% Tomas Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-08-26
function [K,R,T] = m2krt(M)
[K,R,T] = P2KRC(M);