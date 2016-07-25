function c = pcat(p)
%
% p ... a pyramid, where p{i} is a matrix of size k(i) x l(i) x m 
% c ... all in one row in a large matrix for imshow
%
% imshow(pcat(p)) 

for i=1:length(p)                                   % sizes
    s(i,:) = size(p{i});
end
if i==1 s(2,:) = [0 0]; end
e = cumsum(s(2:end,1));
b = [1;1;e+1];
b = b(1:end-1);
e = [s(1,1);e];

c = zeros(s(1,1),sum(s(1:2,2)),size(p{1},3));   % composition image
c(b(1):e(1),1:s(1,2),:) = p{1};                 % compose
for i=2:length(p)                               
    c(b(i):e(i),s(1,2)+[1:s(i,2)],:) = p{i};
end
return