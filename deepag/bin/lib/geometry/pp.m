function s = pp(x)
% 
% x ...  3 x n matrix [a b c; ...]'
% s ... { x(:,i)*x(:,i)' ... }

if size(x,2)==1
    s{1} = x*x';   
else
    for i=1:size(x,2)
        s{i} = x(:,i)*x(:,i)';   
    end
end