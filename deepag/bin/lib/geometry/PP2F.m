% F = PP2F(P1,P2) - Fundamental matrix F from projections matrices P1, P2
%
% P1, P2 = projection matrices
% F      = Fundamental matrix F = [P2*null(P1)]_x * P2 * pinv(P1) scaled such that ||F||_F = sqrt(2)*||C2-C1||

% T.Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-07-11
function F = PP2F(P1,P2)
if all(isfinite(P1(:))) && all(isfinite(P2(:)))
    % camera centers
    C1 = h2a(null(P1));
    C2 = h2a(null(P2));
    % fundamental matrix
    F = xx(P2*null(P1))*P2*pinv(P1);
    F = F/vnorm(F(:))*vnorm(C2-C1)*sqrt(2);
else 
    F = nan(3); % to provide uniform output in cellrrays of P's where some P's are not defined
end
