% hms = sec2hms(s) - seconds s to 'h/m/s' string hms

% 2015-05-10 pajdla@cmp.felk.cvut.cz
function hms = sec2hms(s)

s = round(s);
h = fix(s/3600); s = rem(s,3600);
m = fix(s/60); s = rem(s,60);

hms = sprintf('%02d h %02d m %02d s',h,m,s);  