% v = trandn(sz,sigma,min,max[,nit]) - truncated randn generator such that min <= c <= max
%
% sz       = array of sizes passed to randn(size)
% sigma    = sigma used as sigma*randn(size)
% min, max = truncation limits  
% nit      = number of itterations
% v        = random number matrix of size sz

% (c) T.Pajdla, www.neovision.cz, Sep 22, 2004
function v = trandn(sz,sigma,mn,mx,nit)

if nargin<5
    nit = 1000;
end

v  = sigma*randn(sz);  
mi = v<mn | v>mx;             % elements out of bounds
for i=1:nit
    if ~any(mi(:)) break; end % all in bounds     
    vs    = sigma*randn(sz);  % generate new sample
    v(mi) = vs(mi);           % get new values where out of bounds  
    mi    = v<mn | v>mx;      % elements out of bounds
end
if i==nit
    error(sprintf('trand - %d iterations without convergence, %d elements out of bounds',nit,sum(mi(:))));
end

