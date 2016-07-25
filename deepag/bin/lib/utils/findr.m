% ix = findr(A,r) - find row r in a matrix (cellarray od vectors) A
%
% A = m x n matrix
% r = 1 x n vector
% ix = index such that A(ix,:)==ones(size(ix,1),1)*r

% 2015-09-06 pajdla@cmp.felk.cvut.cz
function ix = findr(A,r)

if iscell(A)
    r = r(:);
    ix = cellfun(@(x) all(x(:)==r),A);
else
    r = r(:)';
    ix = all((A-(ones(size(A,1),1)*r))==0,2);
end
ix = find(ix);