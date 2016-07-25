function updtcrc(obj,eventdata,h) 

u = get(gcf,'CurrentPoint');

x  = get(h,'xdata');
y  = get(h,'ydata');
xc = mean(x);
yc = mean(y);
x  = x - xc;
y  = y - yc;
r  = sqrt(x.^2+y.^2);
r  = mean(r);
if r>eps
    x  = x/r;
    y  = y/r;
end
rr = sqrt((xc-u(1))^2+(yc-u(2))^2);
x  = rr*x+xc; 
y  = rr*y+yc; 
set(h,'xdata',x,'ydata',y);
