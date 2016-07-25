% ix = valfind(v,x) find posiiton index of occurences of elements of v in x
% v  ... a vector of n values
% x  ... a vector of m values
% ix ... {[i,j],{k},} n cells, ix{i} are positions of v(1) in x

% Tomas Pajdla (pajdla@ciirc.cvut.cz)
% 2013-11-01
function ix = valfind(v,x)

ix = cell(1,numel(v));
for i=1:numel(v)
    ix{i} = find(v(i)==x);
end
    
    
