% [B,d] = genp(A[,tol]) - Gaussian elimination - no pivoting
%
% A   ... matrix
% tol ... tollerance, max(size(A))*eps(class(A))*norm(A,'inf') implicit
% B   ... eliminated matrix
% d   ... removed rows
%
% See also RREF

% T. Pajdla, 13 June 2006
% pajdla@neovision.cz

function [B,d] = genp(A,tol)
[m,n] = size(A);

% Compute the default tolerance if none was provided.
if (nargin < 2), tol = max(size(A))*eps(class(A))*norm(A,'inf'); end
if abs(A(1,1))<tol 
    B = [];d = [];
    warning('genp: A(1,1) = 0');
    return
else
    m = 0; d = []; % no linearly dependent rows yet
    k = 1;
    h(1) = 1; % heads = the first non-zero elements in row k
    B(1,:) = A(1,:); % the first row must be there
end
for i=2:size(A,1) % for all remaining rows of A
    b = A(i,:);
    for j=1:min(k,min(size(A))) % eliminate by rows in B
        b = b - b(h(j))*(B(j,:)/B(j,h(j)));
    end
    % find the first non-zero element left
    nix = find(abs(b)>tol); 
    if isempty(nix) % b completely eliminated -> dependent -> remember row index
        m = m + 1;
        d(m) = i;
    else % not dependent -> reduce & add to B on the proper place
        nix = nix(1);
        b(1:nix-1) = 0; % zero all what is below tol
        b = b/norm(b,2);
        k = k + 1;
        B(k,:) = b;
        h(k) = nix; 
        [hs,hi] = sort(h); % make B triangular again
        B = B(hi,:);
        h = hs;
    end
end

