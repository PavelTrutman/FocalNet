% c = colorlist(t,i) - colors for plotting
%
% t = type of output
%     'txt' ... c ={'r','g',...} (missing)
%     'rgb' ... c ={[0 0 1],...}
%     'rnd' ... c ={[. . .],...} of 64 same randomly chosen colors
% i = index color modulo the color size
%     number ... returns color
%     vector ... returns cell array of colors

% pajdla@cmp.felk.cvut.cz
% 2015-05-17
function c = colorlist(t,i)
if nargin<1
    t = 'txt';
end
switch t
    case 'txt'
        c = {'b','g','r','c','m','w','k','y'};
    case 'rgb'     
        c = {[0.00 0.00 1.00]
             [0.00 1.00 0.00]
             [1.00 0.00 0.00]
             [0.00 1.00 1.00]
             [1.00 1.00 0.00]
             [1.00 0.00 1.00]
             [1.00 1.00 1.00]
             [0.00 0.00 0.00]};
    case 'rnd'
        c = lines;
        s = rng; rng(pi);
        c(8:end,:) = rand(size(c,1)-8+1,3); 
        c = mat2cell(c,ones(size(c,1),1),3);
        rng(s);
    otherwise
        error('%s - unknown value',t);
end
% select the colors
if nargin<2
    i = 1:numel(c);
end
i = mod(i-1,numel(c))+1;
if numel(i)>1
    c = c(i);
else
    c = c{i};
end
