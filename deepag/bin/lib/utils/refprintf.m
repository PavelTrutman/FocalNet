% nbytes = refprintf(vargin) - rewrite fprintf over
%
% varargin, nbytes ... nbytes = fprintf(varargin)

% 2015-05-10 pajdla@cmp.felk.cvut.cz
function nbytes = refprintf(varargin)

switch nargin
    case 1
        fprintf(repmat('\b',1,length(varargin{1})));
        if nargout>0
            nbytes = fprintf(varargin{1});
        else
            fprintf(varargin{1});
        end
    otherwise
        s = sprintf(varargin{2:end});
        fprintf(repmat('\b',1,length(s)));
        if nargout>0
            nbytes = fprintf(varargin{1},s);
        else
            fprintf(varargin{1},s);
        end
end
