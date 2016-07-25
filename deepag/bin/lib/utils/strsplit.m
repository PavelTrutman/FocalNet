% sp = strsplit(s,d) - split a string into substrings delimited by d
% 
% s  ... a string
% d  ... a delimitor see strtok
% sp ... a cell array of substrings
function sp = strsplit(s,d)

for i=1:length(s)
   [si, s] = strtok(s,d);
   if isempty(si), break;  end
   sp{i} = si;
end