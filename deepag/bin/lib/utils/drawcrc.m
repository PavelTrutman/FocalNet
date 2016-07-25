function h = drawcrc(p) 

fi = 0:1/20:2*pi; 
x  = cos(fi); 
y  = sin(fi);

h = plot(p(3)*x+p(1),p(3)*y+p(2),'m');
