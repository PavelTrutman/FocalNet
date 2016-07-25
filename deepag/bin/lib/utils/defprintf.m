% nbytes = refprintf(vargin) - delete space for fprintf 
%
% varargin, nbytes ... nbytes = fprintf(varargin)

% 2015-08-23 pajdla@cmp.felk.cvut.cz
function nbytes = defprintf(varargin)

switch nargin
    case 1
        nbytes = length(varargin{1});
        fprintf(repmat('\b',1,nbytes));
    otherwise
        s = sprintf(varargin{2:end});
        nbytes = length(s);
        fprintf(repmat('\b',1,nbytes));
end
