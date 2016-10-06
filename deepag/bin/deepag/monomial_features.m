% calculate feature vector with coefficients corresponding to monomials in
% some equation, see ShareLatex
%
% (Zuzana?)
function [Mvec] = monomial_features(x, xp)
  M = [-xp(:,2), x(:,2), xp(:,1), x(:,1), -x(:,1).*xp(:,2) - x(:,2).*xp(:,1), 2*x(:,1).*xp(:,1) - 2*x(:,2).*xp(:,2), 2*x(:,1).*xp(:,2), -2*x(:,1), 2*x(:,2), 2*ones(size(x,1),1), -2*x(:,1).*xp(:,1), -2*xp(:,2), -2*ones(size(x,1),1), -2*x(:,2).*xp(:,2), 2*x(:,2).*xp(:,1), 2*xp(:,1), -xp(:,2), -x(:,2), xp(:,1), -x(:,1), x(:,1).*xp(:,2) + x(:,2).*xp(:,1), -2*ones(size(x,1),1), 2*x(:,2).*xp(:,2), -2*x(:,2).*xp(:,1), 2*xp(:,1), 2*x(:,1).*xp(:,2), -2*ones(size(x,1),1), -2*x(:,1).*xp(:,1), 2*xp(:,2), xp(:,2), x(:,2), -xp(:,1), -x(:,1), x(:,1).*xp(:,2) - x(:,2).*xp(:,1), 2*x(:,1), 2*x(:,2), - 2*x(:,1).*xp(:,1) - 2*x(:,2).*xp(:,2), xp(:,2), -x(:,2), -xp(:,1), x(:,1), x(:,2).*xp(:,1) - x(:,1).*xp(:,2)];
  Mr = rref(M);
  if rank(Mr(:, 1:7)) < 7
    Mvec = NaN;
    return;
  end
  Mrs = Mr(:, 8:end);
  %Mvec = Mrs([4 9 15:29 36:56 59 64 72 80 88 96 99:119 122 127:148 155 163 171 179 183:189 193 198 204:211 219 227 235 239:245]);
  Mvec = Mrs(:);
end