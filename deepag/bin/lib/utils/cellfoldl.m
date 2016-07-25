% o = cellfoldl(@f,c) - f(c{1}, ... f(c{end-2},f(c{end-1},c{end})))
%
% @f = reference to a function f
% c  = cell with at least 2 elements
% o  = f(c{1}, ... f(c{end-2},f(c{end-1},c{end}))) (by reccurent calls)

% pajdla@cmp.felk.cvut.cz, 2015-05-08
function o = cellfoldl(f,c)

if numel(c)==2
    o = f(c{1},c{2});
elseif numel(c)>2
    o = f(c{1},cellfoldl(f,c(2:end)));
else
    error('c must have at least two elements!');
end
