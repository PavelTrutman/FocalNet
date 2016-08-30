% [l,r] = l2Dfit(y) - Line fit
% y ... 2 x n points
% l ... 3 x 1 line coordinates
% r ... 1 x n residuals
% z ... 4 x n [the closest points z on the line; y]
% See l2Dres

% (c) T.Pajdla, www.cvut.cz, 2016-08-08
function  [l,r,z,C] = l2Dfit(y,~,dbg)

if nargin>0
    if nargin<3
        dbg = false;
    end
    % my = mean(y,2); slow
    n = size(y,2);
    my = sum(y,2)/n;
    if dbg % straghtforward formulation but slower
        z  = y - my*ones(1,n);
        C  = z*z';
    else % faster
        % C = C-n*(my*my'); % (my*my') - parenthesis needed to make sure that it is symmetric - still slower than the formula below
        C = y*y'-n*[my(1)*my(1) my(1)*my(2); my(1)*my(2) my(2)*my(2)];
    end
    [v,e] = eig(C);
    l  = [v(:,1); -v(:,1)'*my];
    if nargout>1
        r = l2Dres(l,y);
    end
    if nargout>2
        z = [y-l(1:2)*(l'*a2h(y));y];
    end
else % unit tests
    % Test 1
    x = rand(2,1000);
    [l1,~,~,C1] = l2Dfit(x);
    [l2,~,~,C2] = l2Dfit(x,[],true);
    l(1) = vnorm(l2-l1)<1e-6;
    % Test 2
    x = [1 2 3;2 3 4];
    l1 = l2Dfit(x);
    l(2) = all(abs(l1'*a2h(x))<1e-10); 
    l(3) = vnorm(l1-[-1;1;-1]/sqrt(2))<1e-6;
end