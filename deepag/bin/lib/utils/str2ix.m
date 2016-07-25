% function i = str2ix(c,s) - index of string occurences in a celarray
% c = celarray of strings
% s = celarray of strings
% i = celarray of indices of occurences of elements of s in c

function i = str2ix(c,s)
for j=1:length(s)
  i{j} = find(~cellfun(@isempty,strfind(c,s{j})));
end