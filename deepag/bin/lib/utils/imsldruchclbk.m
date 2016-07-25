% imsldruchclbk() - image slider callback
%
% See imslider

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 3 Nov 2009
function imsldruchclbk(a,b,c)

iminf = get(gcf,'userdata');
n = iminf.n;
nb=floor(n/2);
nt=ceil(n/2);
p = round(get(iminf.uch,'Value'));
b = p-nb; t = p+nt;
if t>size(iminf.lpos,2) 
    t=size(iminf.lpos,2); 
    b=t-n; 
end
if b<1
    t = t + (1-b);
    b = b + (1-b);
    t = min(t,size(iminf.lpos,2));
end
axis([iminf.lpos([b t]) 1 iminf.imh]);
return