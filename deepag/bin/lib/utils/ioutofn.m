% s = ioutofn(i,n[,maxn]) - fprintf '[i/n]'  
%
% i = current value
% n = out of n
% maxn = maximal value of n to determine the width of printing (= n if missing)
% s = string

% 2015-05-10 pajdla@cmp.felk.cvut.cz
function s = ioutofn(i,n,maxn)
if nargin<3 maxn = n; end
N = length(sprintf('%d',maxn)); % max number of digits
s = sprintf(sprintf('[%%0%dd/%%0%dd]',N,N),i,n);

