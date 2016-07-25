-- Van der Waerden Syzygy 
-- T. Pajdla (pajdla@cvut.cz)
-- B. Sturmfels. Algorithms in Invariant Theory. Springer 2008.

-- Constraints on 3x3 minors of 3x6 matrix
-- R = ZZ/31991[p_(1,1)..p_(3,6),d_{1, 2, 3}, d_{1, 2, 4}, d_{1, 3, 4}, d_{2, 3, 4}, d_{1, 2, 5}, d_{1, 3, 5}, d_{2, 3, 5}, d_{1, 4, 5}, d_{2, 4, 5}, d_{3, 4, 5}, d_{1, 2, 6}, d_{1, 3, 6}, d_{2, 3, 6}, d_{1, 4, 6}, d_{2, 4, 6}, d_{3, 4, 6}, d_{1, 5, 6}, d_{2, 5, 6}, d_{3, 5, 6}, d_{4, 5, 6},MonomialOrder=>Lex]
clearAll
m = 2;
n = 8;
o = apply(m,i->1);
ix = subsets(1..n,m);
R = ZZ/31991[p_(1,1)..p_(m,n),(i->d_i) \ ix,MonomialOrder=>Lex];
M = transpose genericMatrix(R,p_(1,1),n,m);
P = (i->d_i-det(M_(i-o))) \ ix;
I = ideal P;
numgens R, #P, dim I, degree I
J = gens gb I;

-- Move to a new ring with d_{} variable only and define there an ideal
-- given yb relations between original minors
Rd = ZZ/31991[(i->d_i) \ ix,MonomialOrder=>Lex];
Pd = matrix J_{0..153}; -- J_{0..59}; -- J_{0..5}; -- J_{0..671}
Id = sub(ideal Pd,Rd);
numgens Rd, #Pd, dim Id, degree Id

A = toList(1..numgens Rd);
L = (i->d_(ix_i)-A_i) \ {0,1,2,3,4,6,7,10,11,15,16,21,22}
IdL = Id + ideal L;
B = gens gb IdL

