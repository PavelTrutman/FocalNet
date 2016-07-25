% f = plotq(q[,opts]) - plot conic q
% 
% q    ... 3x3 conic
% f    ... figure handle returned by plot
% opts ... = []        ... f = point coordinates
%            otherwise ... plot options plot(...,...,opts)

% (c) T.Pajdla, www.neovision.cz, Oct 25, 2004
function f = plotq(q,opts)

if nargin < 2
    opts = ' ';
end

x  = sampleq(q);
x  = h2a(x);
ix = all(isfinite(x));
x  = x(:,ix);
if nargin<2
        f = plot(x(1,[1:end 1]),x(2,[1:end 1]));
else
    if ~isempty(opts)
        f = plot(x(1,[1:end 1]),x(2,[1:end 1]),opts);
    else
        f = x;
    end

end
