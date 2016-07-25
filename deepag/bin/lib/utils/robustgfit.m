% [m,s,d,j] = robustgfit(t,x,n,mi,si,tol,sim) - Robust Guassian distribution fit
%
% t     ... samples
% x     ... abscisa
% n     ... frequencies
% mi    ... initial mean
% si    ... imitial std
% tol   ... tolerance (to multiply s)
% sim   ... {'Likelihood','DistrDiff'}
% m     ... mean
% s     ... sigma
% d     ... similarity of the distribution and the data
% j     ... iterations

% T. Pajdla, pajdla@neovision.cz
% 30 Dec 2006
function [m,s,d,j] = robustgfit(t,x,n,mi,si,tol,sim)

N = 10; % max num of interations
m = zeros(1,N);
s = zeros(1,N);
% initial Gaussian model
m(1) = mi;
s(1) = si;
switch sim % agreement on inliers
    case 'Likelihood'        
        d(1:2,1) = sum((smplgaus(t,m(1),s(1))))/numel(t); % agreement on inliers
    case 'DistrDiff'
        ixx = abs((x-m(1)))<tol*s(1);
        d(1:2,1) = sum(-abs((n/sum(n)-smplgaus(x,m(1),s(1))).*ixx));
end
for j=2:N
    % new distribution
    ix = abs((t-m(j-1)))<tol*s(j-1);
    m(j) = mean(t(ix));
    s(j) = std(t(ix));
    % compare the quality of distributions
    switch sim % agreement on inliers
        case 'Likelihood'
            d(1:2,j) = [sum((smplgaus(t(ix),m(j-1),s(j-1))))
                        sum((smplgaus(t(ix),m(j),s(j))))]/numel(t);
        case 'DistrDiff'
            ixx = abs((x-m(j-1)))<tol*s(j-1);
            d(1:2,j) = [sum(-abs((n/sum(n)-smplgaus(x,m(j-1),s(j-1))).*ixx))
                        sum(-abs((n/sum(n)-smplgaus(x,m(j),s(j))).*ixx))];
    end
    if d(2,j) <= d(1,j)
        j = j-1;
        break;
    end
end
m = m(j);
s = s(j);
d = d(2,j);


