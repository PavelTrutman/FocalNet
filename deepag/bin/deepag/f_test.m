f = zeros(1, size(pdAll.C.cam, 2));

for i = pdAll.C.cix
  f(i) = pdAll.C.cam{i}.f;
end
inliers = pdAll.C.cix;

inliersLast = size(inliers, 2) + 1;
while (inliersLast - size(inliers, 2)) ~= 0
  inliersLast = size(inliers, 2);
  outliers = find(f > (mean(f(inliers)) + 3*std(f(inliers))));
  fprintf(['Removed ', num2str(size(outliers, 2)), 'cams.\n']);
  inliers = setdiff(inliers, outliers);
end