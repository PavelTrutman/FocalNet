% {P1,P2} = E2PP(E) - Essential Matrix E to projection matrices P1, P2 (4 solutions)
%
% E = essential matrix
% P1 = projection matrix = {[I 0],[I 0],[I 0],[I 0]}
% P2 = projection matrix = {     ,     ,     ,     }
%
% such that
%
%  x2^T E x1 = 0 => a1 x1 = P1 * X & a2 x2 = P2 X
%
% The distance between the camera centers is set to ||E||_F/sqrt(2).

% T. Pajdla, pajdla@cmp.felk.cvut.cz, 2015-07-11
% http://cmp.felk.cvut.cz/~pajdla/gvg/GVG-2015-Lecture.pdf (12.4.1 Camera Computation)
function P12 = E2PP(E)
if nargin>0
    % get size of t, which is encoded as ||E||_F = sqrt(2)*||t||
    nt = vnorm(E(:))/sqrt(2);
    % normalize E
    G = E/nt; % 12.53
    % get the translation direction
    v = null3x3r2(G); %12.54
    v = v/norm(v);
    % rotation matrices
    V = xx(v); % 12.57
    g1 = G(:,1); g2 = G(:,2); g3 = G(:,3);
    v1 = V(:,1); v2 = V(:,2); v3 = V(:,3);
    % s = +1
    B = [g1 g2 g3 cr(g1,g2) cr(g2,g3) cr(g1,g3)];
    A = [v1 v2 v3 cr(v1,v2) cr(v2,v3) cr(v1,v3)];
    R{1} = B*pinv(A);
    % R{1} = XY2rot(A,B); % mumerically exact rotation -> this removes some good matches - why?
    % s = -1
    B = B*diag([-1 -1 -1 1 1 1]);
    R{2} = B*pinv(A);
    % R{2} = XY2rot(A,B); % mumerically exact rotation -> this removes some good matches - why?
    % translation vectors
    v = nt*v; % scale the translation
    t{1} = +v;
    t{2} = -v;
    % compose P's
    I = [1 0 0;0 1 0;0 0 1];
    O = [0;0;0];
    P1 = {[I O],[I O],[I O],[I O]};
    P2 = {R{1}*[I -t{1}],R{1}*[I -t{2}],R{2}*[I -t{1}],R{2}*[I -t{2}]};
    P12 = {P1,P2};
else % unit tests
    test.test = {'coordinate system orientation'};
    P1 = eye(3,4);
    P2 = [eye(3) -[1;0;0]];
    E  = PP2F(P1,P2);
    PP = E2PP(E);
    d = cellfun(@(p) max(abs(p(:)-P2(:))),PP{2});
    test.ok = any(d<1e-8);
    P12 = test;
end

function z = cr(x,y)
  z = [x(2)*y(3)-x(3)*y(2);-x(1)*y(3)+x(3)*y(1);x(1)*y(2)-x(2)*y(1)];
return