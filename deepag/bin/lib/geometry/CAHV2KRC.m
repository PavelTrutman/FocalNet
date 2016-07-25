% [K,R,C] = CAHV2KRC(c) - CAHV to KRC camera model conversion
%
% c       = CAHV camera model 
%           c.C, c.A, c.H, c.V, c.b, c.type = 'CAHV', see X2u.m
% K, R, C = KRC camera matrices, see X2u.m
 

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-07-12
function [K,R,C] = CAHV2KRC(c)

hc = c.A'*c.H;
vc = c.A'*c.V;
hs = vnorm(xx(c.A)*c.H);
vs = vnorm(xx(c.A)*c.V);
Hp = (c.H - hc*c.A)/hs; 
Vp = (c.V - vc*c.A)/vs;
R  = [Hp';Vp';c.A'];
C  = c.C;
K  = [hs 0 hc;0 vs vc;0 0 1];
K  = K/K(2,2)/c.b(2); % get 1/bx 1/by 1/f in m on the diagonal of K
