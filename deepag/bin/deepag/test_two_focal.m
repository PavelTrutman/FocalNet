d = zeros(245, 3*10000*100);
l = 1;
for i = 1:3
for j = randperm(size(dataset{i}.cameraPairs, 2))
  for k = randperm(size(dataset{i}.cameraPairs{j}.Xmask, 2))
    d(:, l) = two_focal(pdAll.C.cam{dataset{i}.cameraPairs{j}.cix1}.u(:, dataset{i}.cameraPairs{j}.u1mask{k})', pdAll.C.cam{dataset{i}.cameraPairs{j}.cix2}.u(:, dataset{i}.cameraPairs{j}.u2mask{k})');
    l = l + 1;
  end
end
end

nze = (sum(d ~= 0, 2) ~= 0);

clear i j k l;