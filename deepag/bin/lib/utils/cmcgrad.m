% g = cmcgrad(imx,diffMask)
% 
% Gradient for camcheck
%
% See CAMCHECK

% T. Pajdla, pajdla@neovision.cz
% 30 Dec 2006

function g = cmcgrad(imx,diffMask)

% border size
s = floor(length(diffMask)/2); 
% gradient magnitude
g = sqrt(conv2(imx,diffMask,'same').^2 + conv2(imx,diffMask','same').^2)/(sum(abs(diffMask))/2);
% fill borders with repeated valid values
g(1:s,:) = repmat(g(s+1,:),[s 1]); % top
g(end-(0:s-1),:) = repmat(g(end-s,:),[s 1]); % botom
g(:,1:s) = repmat(g(:,s+1),[1 s]); % left
g(:,end-(0:s-1)) = repmat(g(:,end-s),[1 s]); % left