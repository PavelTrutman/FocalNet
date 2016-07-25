%url = fullurl(varargin) - fullurl consruction  (see fullfile)
function u = fullurl(varargin)
u = fullfile(varargin{:});
u = strrep(u,'\','/');