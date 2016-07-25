% h = plotCam(P or K,R,T) - plot camera
% 
% P or K,R,T ... 3x4 or 3x3 camera matrices
% h          ... plot handles

% (c) T.Pajdla, www.neovision.cz, Sept 14, 2007
function h = plotCam(K,R,T)

if ~isempty(K)
    if nargin>2
        K = [K*R -K*R*T];
    end
    h = camplot(K);
else
    h = [];
end

