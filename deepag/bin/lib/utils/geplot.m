% geplot(A[,X,L,opt])- plot a graph 
%
% A = n x n edge weight matrix, i.e. A(i,j) = weight of edge i-j
% X = 2(3) x n vertex coordinates, regular n-gon if =[] or missing
% L = labels matrix of strings or a vector of numbers, can be omited
% opt = string of vertex plot options

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 27 Oct 2009
function geplot(A,X,L,opt) 

if size(A,1)~=size(A,2)
    error('A must be square');
end
n = size(A,1);
if nargin<4
    opt='r.';
end
if nargin<3
    L = [];
end
if nargin<2 || isempty(X)
    t = 0:2*pi/n:2*pi-2*pi/n;
    X = [cos(t) ; sin(t)]; 
end

ix = sub2ind(size(A),1:size(A,1),1:size(A,2));
A(ix) = 0;
A = triu(A);
ix = find(A);
[i,j]=ind2sub(size(A),ix);
switch size(X,1)
    case 2 
        plot(X(1,:),X(2,:),opt);
    case 3 
        plot3(X(1,:),X(2,:),X(3,:),opt);
    otherwise
        error('X must be 2 or 3 time n');
end
if ~isempty(i)
    set(gca,'color',[0,0,0]);
    %
    c = A(ix);
    c = c-min(c);
    c = c/max(c);
    ix = find(~isfinite(c)); % unknown similarity
    c = [1;1;1]*c';
    if ~isempty(ix)
        c(:,ix) =  [1;0;0]*ones(size(ix)); % unknown similarity mark red
    end
    c0 = 0.3;
    switch size(X,1)
        case 2
            for k=1:length(i)
                line([X(1,i(k));X(1,j(k))],[X(2,i(k));X(2,j(k))],'color',(1-c0)*c(:,k)+c0);
            end
        case 3 
            for k=1:length(i)
                line([X(1,i(k));X(1,j(k))],[X(2,i(k));X(2,j(k))],[X(3,i(k));X(3,j(k))],'color',(1-c0)*c(:,k)+c0);
            end
    end    
end 
if ~isempty(L)
    if ~isstr(L)
        L=num2str(L(:));
    end
    hs = ishold;
    if ~hs, hold on, end
    switch size(X,1)
        case 2
            text(X(1,:),X(2,:),L);
        case 3
            text(X(1,:),X(2,:),X(3,:),L);
    end
    if ~hs, hold off, end
end
end
