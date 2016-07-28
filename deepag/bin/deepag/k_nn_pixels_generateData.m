
valSize = 500;
overMarginCoef = 1.05;

v = [repmat([-6/7; -4/7; -2/7; 0; 2/7; 4/7; 6/7], 2, 1); -repmat([-6/7; -4/7; -2/7; 0; 2/7; 4/7; 6/7], 2, 1)];

u_tr = v;
u_val = repmat(v, 1, valSize) + unifrnd(-1/7*overMarginCoef, 1/7*overMarginCoef, 28, valSize);

save('pixels_syntetic.mat', 'u_tr', 'u_val', '-v7.3');