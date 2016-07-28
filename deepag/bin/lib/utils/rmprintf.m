% S = rmprintf(S) - delete last string in celarray S from stdout
%
% S = cellarray of strings 

% 2015-08-23 pajdla@cmp.felk.cvut.cz
function S = rmprintf(S)

if numel(S)>0
    fprintf(repmat('\b',1,length(sprintf(S{end})))); 
    if numel(S)>1
        S = S(1:end-1);
    else
        S = {};
    end
end