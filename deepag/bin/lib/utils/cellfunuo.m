% o = cellfunuo(varargin) => o = cellfun(varargin,'UniformOutput',false)
%
% @f = reference to a function f
% c  = cell
% o  = cellfun(varargin,'UniformOutput',false)

% pajdla@cvut.cz, 2016-04-22
function o = cellfunuo(varargin)

o = cellfun(varargin{:},'UniformOutput',false);
