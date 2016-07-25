% S = adprintf(S,s) - fprintf s and append to celarray S
%
% S = cellarray of strings 
% s = string

% 2015-08-23 pajdla@cmp.felk.cvut.cz
function S = adprintf(S,s)

fprintf(s);
S{end+1}=s;
drawnow;