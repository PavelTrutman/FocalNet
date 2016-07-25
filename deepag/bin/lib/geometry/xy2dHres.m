% [e,J] = xy2dHres(x,y[,m,opts,X,ix]) - Residuals of homohgraphies & targer points
%
% y{i}  ... 3 x n 2D observed image points
% m     ... {'X','HX'}
%             X ... residuals of given x w.r.t. etimated best H's
%               x ... 2*n x 1 (3 x n) 2D target points
%             HX ... residuals of given x a H's
%               x ... 2*n x 1 + length(y)*9 2D target points & homographies
% opts ... optimization parameters
%          opts = optimset(...)
%          implicit: optimset('Display','none','TolX',1e-5,'MaxIter',100);
% X     ... 2*n x 1 (3 x n) all 2D target points
% ix    ... index of points that were modified (mising & [] = all)
% e     ... residuals
% J     ... Jacobian ... not yet implemented

% (c) T.Pajdla, www.neovision.cz, Feb 1 2006
function [e,J] = xy2dHres(x,y,m,opts,X,ix)
if nargin < 5
    X = x;
end
if nargin<4
    opts  = optimset('Display','none','TolX',1e-5,'MaxIter',100);
end
if nargin < 3
    m = 'X';
end
switch m
    case 'X' % residuals of given x w.r.t. the best H's
        if nargin < 6 | isempty(ix)
            ix = logical(ones(size(x)));
        end
        if size(x,2)==1
            x = reshape(x,2,length(x)/2); % x to column vectors
            x = a2h(x); % homogeneous coordinates 
        end
        X(:,ix) = x; % update points
        for i=1:length(y) % for all measurements
            [H,e{i}]=xy2H(X,y{i},{'DLT','y-y'},0,opts); % best homography fit
        end
        e = cat(1,e{:}); % residuals
        if nargout>1 % compute the Jacobian
            J = [];
            fprintf('\n');
        end
    case 'HX' % residuals of given x a H's
        h = x(end-(length(y)*9-1):end); % get h's
        x = x(1:end-(length(y)*9)); % get x's
        x = reshape(x,2,length(x)/2); % x to column vectors
        x = a2h(x); % homogeneous coordinates 
        if nargin < 6 | isempty(ix)
            ix = logical(ones(1,size(x,2)));
        end
        X(:,ix) = x; % update points
        for i=1:length(y) % for all measurements
            H = reshape(h((i-1)*9+[1:9]),3,3)'; % H stored rowwise 
            % e{i} = xy2Hres(H,X,y{i},'y-y');  
            % e{i} = sqrt(sum((h2a(y{i})-h2a(H*x)).^2))';
            e{i} = reshape(h2a(y{i})-h2a(H*X),2*size(x,2),1); % large scale methods need as many independent variables as possible
        end
        e = cat(1,e{:}); % residuals
        % e = reshape(e,16,length(e)/16);
        % e = sqrt(sum(e.^2));
        if nargout>1 % compute the Jacobian
            % make an empty matrix
            J = sparse([],[],[],2*size(X,2)*length(y),2*size(X,2)+9*length(y),14*size(X,2)*length(y));
            E = spdiags(ones(2*size(X,2),1),0,2*size(X,2),2*size(X,2)); % the identity
            W = sparse(2*size(X,2),9); % empty non-trivial part of J
            for i=1:length(y) % for all measurements
                H = reshape(h((i-1)*9+[1:9]),3,3)'; % H stored rowwise
                w   = H*X;
                w1  = ([1;1;1]*w(1,:));
                w2  = ([1;1;1]*w(2,:));
                w3  = ([1;1;1]*w(3,:));
                yw3 =  X./w3;
                yw13 = w1.*yw3./w3;
                yw23 = w2.*yw3./w3;
                yw3 = yw3';yw13 = yw13';yw23 = yw23';
                % construct the non-trivial part
                W(1:2:end-1,:) = [yw3 sparse(size(yw3,1),size(yw3,2)) yw13]; 
                W(2:2:end,:)   = [sparse(size(yw3,1),size(yw3,2)) yw3 yw23];
                % insert the identity part of J
                J((i-1)*2*size(X,2)+[1:2*size(X,2)],1:2*size(X,2)) = E; 
                % insert the non-trivial part of J
                J((i-1)*2*size(X,2)+[1:2*size(X,2)],2*size(X,2)+(i-1)*9+[1:9]) = W; 
            end
            %% add the identity part
            % surprisingly, the following code took longer
            % i = 1:size(J,1);
            % j = rem(i-1,2*size(X,2))+1;
            % ix = sub2ind(size(J),i,j);
            % J(ix) = 1;
        end        
end
return
