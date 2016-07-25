% i = gesymix(e) - Match symetric edges of a graph
%
% e = edges - n x 2 matrix of vertex numbers
% i = n x 2 index = [i1 i2] such that e(i1,[1 2]) = e(i2,[2 1]), sorted lex i(1)>i(2)

% 2015-08-16 pajdla@cmp.felk.cvut.cz
function ix = gesymix(e)

s = max(e(:))*[1 1]; % size of the adjacency matrix
i1 = sub2ind(s,e(:,1),e(:,2)); % index of e
i2 = sub2ind(s,e(:,2),e(:,1)); % index of tansposed e
[si1,six1]=sort(i1);
[si2,six2]=sort(i2);
if all(si1==si2)
    ix = [six1 six2];
    ix = sortrows(ix')';
else
    error('gesymix: e is not symmetric')
end