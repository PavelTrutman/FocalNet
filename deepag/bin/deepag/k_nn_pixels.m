function k_nn_pixels()
  % k-nearest neighbours in the normalized image coordinates
  % 
  % Pavel Trutman
  % INRIA, 2016

  deepagpaths;

  % prepare data
  load('../../data/paris/correspondences.mat');
  
  % properties
  margin = 0.2;
  
  % find cams that their correspondences are far from each other
  fprintf('Choosing cameras with correspondences that are far from each other.\n'); tic;
  repS = adprintf({}, 'Training dataset: ');
  cix_tr = findGoodCams(corr_tr.u, margin);
  rmprintf(repS);
  repS = adprintf({}, 'Validating dataset: ');
  cix_val = findGoodCams(corr_val.u, margin);
  rmprintf(repS);
  fprintf('\b\b ... done %s.\n', sec2hms(toc));

  u_tr = corr_tr.u(:, cix_tr);
  %f_tr = corr_tr.f(:, cix_tr);
  %n_tr = corr_tr.norm(:, cix_tr);
  u_val = corr_val.u(:, cix_val);
  %f_val = corr_val.f(:, cix_val);
  %n_val = corr_val.norm(:, cix_val);
  
  % nearest neighbour for points validating dataset
  fprintf('Searching for nearest neighbour for points validating dataset.\n'); tic;
  
  u1x_tr = u_tr(1:7, :);
  u1y_tr = u_tr(8:14, :);
  u2x_tr = u_tr(15:21, :);
  u2y_tr = u_tr(22:28, :);
  
  for i = 1:size(u_val, 2)
    fprintf([num2str(i), '/', num2str(size(u_val, 2)), '\n']);
    u1x_val = u_val(1:7, i);
    u1y_val = u_val(8:14, i);
    u2x_val = u_val(15:21, i);
    u2y_val = u_val(22:28, i);
    
    %dist = Inf(size(u_tr, 2), 1);
    
    % not swapped
    order1noswap = zeros(7, size(u_tr, 2));
    order2noswap = zeros(7, size(u_tr, 2));
    noswap = true(1, size(u_tr, 2));
    for j = 1:7
      hxBound = repmat(u1x_val(j, :) + margin, 7, sum(noswap));
      hyBound = repmat(u1y_val(j, :) + margin, 7, sum(noswap));
      lxBound = repmat(u1x_val(j, :) - margin, 7, sum(noswap));
      lyBound = repmat(u1y_val(j, :) - margin, 7, sum(noswap));
      match = ((u1x_tr(:, noswap) < hxBound) & (u1x_tr(:, noswap) > lxBound) & (u1y_tr(:, noswap) < hyBound) & (u1y_tr(:, noswap) > lyBound));
      tmp = zeros(size(noswap));
      tmp(noswap) = (sum(match, 1) == 1);
      order1noswap(j, noswap) = findFirstInCol(match);
      noswap = noswap & tmp;
    end
    for j = 1:7
      hxBound = repmat(u2x_val(j, :) + margin, 7, sum(noswap));
      hyBound = repmat(u2y_val(j, :) + margin, 7, sum(noswap));
      lxBound = repmat(u2x_val(j, :) - margin, 7, sum(noswap));
      lyBound = repmat(u2y_val(j, :) - margin, 7, sum(noswap));
      match = ((u2x_tr(:, noswap) < hxBound) & (u2x_tr(:, noswap) > lxBound) & (u2y_tr(:, noswap) < hyBound) & (u2y_tr(:, noswap) > lyBound));
      tmp = zeros(size(noswap));
      tmp(noswap) = (sum(match, 1) == 1);
      order2noswap(j, noswap) = findFirstInCol(match);
      noswap = noswap & tmp;
    end
    %check uniqueness
    tmp = zeros(size(noswap));
    tmp(noswap) = all(sort(order1noswap(:, noswap)) == repmat((1:7)', 1, sum(noswap)));
    noswap = noswap & tmp;
    tmp = zeros(size(noswap));
    tmp(noswap) = all(sort(order2noswap(:, noswap)) == repmat((1:7)', 1, sum(noswap)));
    noswap = noswap & tmp;
    fprintf([num2str(sum(noswap)), '\n']);

    %swapped
    order1swap = zeros(7, size(u_tr, 2));
    order2swap = zeros(7, size(u_tr, 2));
    swap = true(1, size(u_tr, 2));
    for j = 1:7
      hxBound = repmat(u1x_val(j, :) + margin, 7, sum(swap));
      hyBound = repmat(u1y_val(j, :) + margin, 7, sum(swap));
      lxBound = repmat(u1x_val(j, :) - margin, 7, sum(swap));
      lyBound = repmat(u1y_val(j, :) - margin, 7, sum(swap));
      match = ((u2x_tr(:, swap) < hxBound) & (u2x_tr(:, swap) > lxBound) & (u2y_tr(:, swap) < hyBound) & (u2y_tr(:, swap) > lyBound));
      tmp = zeros(size(swap));
      tmp(swap) = (sum(match, 1) == 1);
      order1swap(j, swap) = findFirstInCol(match);
      swap = swap & tmp;
    end
    for j = 1:7
      hxBound = repmat(u2x_val(j, :) + margin, 7, sum(swap));
      hyBound = repmat(u2y_val(j, :) + margin, 7, sum(swap));
      lxBound = repmat(u2x_val(j, :) - margin, 7, sum(swap));
      lyBound = repmat(u2y_val(j, :) - margin, 7, sum(swap));
      match = ((u1x_tr(:, swap) < hxBound) & (u1x_tr(:, swap) > lxBound) & (u1y_tr(:, swap) < hyBound) & (u1y_tr(:, swap) > lyBound));
      tmp = zeros(size(swap));
      tmp(swap) = (sum(match, 1) == 1);
      order2swap(j, swap) = findFirstInCol(match);
      swap = swap & tmp;
    end
    %check uniqueness
    tmp = zeros(size(swap));
    tmp(swap) = all(sort(order1swap(:, swap)) == repmat((1:7)', 1, sum(swap)));
    swap = swap & tmp;
    tmp = zeros(size(swap));
    tmp(swap) = all(sort(order2swap(:, swap)) == repmat((1:7)', 1, sum(swap)));
    swap = swap & tmp;
    fprintf([num2str(sum(swap)), '\n']);
      
    % points are close
    %{
    if swap
      [u2_tr, u1_tr] = deal(u1_tr, u2_tr);
    end
        
    dist(k) = max([sum((u1_tr(:, order1) - u1_val).^2), sum((u2_tr(:, order2) - u2_val).^2)]);
      
    [m, idx] = min(dist);
    fprintf([num2str(i), '/', num2str(size(u_val, 2)), ': ', num2str(m), '\n']);
    %}
    
  end
  
  
end

%% Functions
  
function cix = findGoodCams(u, margin)
  cix = false(1, size(u, 2));
  
  repS = {};
  for i = 1:size(u, 2)
    if mod(i, 1000) == 0
      repS = rmprintf(repS);
      repS = adprintf(repS, [num2str(i), '/', num2str(size(u, 2))]);
    end
    u1 = [u(1:7, i)'; u(8:14, i)'];
    u2 = [u(15:21, i)'; u(22:28, i)'];
    
    skip = false;
    for j = 1:6
      hBound = repmat(u1(:, j) + repmat(margin, 2, 1), 1, 7 - j);
      lBound = repmat(u1(:, j) - repmat(margin, 2, 1), 1, 7 - j);
      if any(all((u1(:, (j+1):7) < hBound) & (u1(:, (j+1):7) > lBound)))
        skip = true;
        break;
      end
    end
    if skip
      continue;
    end
    for j = 1:6
      hBound = repmat(u2(:, j) + repmat(margin, 2, 1), 1, 7 - j);
      lBound = repmat(u2(:, j) - repmat(margin, 2, 1), 1, 7 - j);
      if any(all((u2(:, (j+1):7) < hBound) & (u2(:, (j+1):7) > lBound)))
        skip = true;
        break;
      end
    end
    if skip
      continue;
    end
    cix(i) = true;
    
  end
  
  rmprintf(repS);
    
end

function b = findFirstInCol(a)
  a = cumsum(a);
  [row, col] = find(a == 0);
  b = zeros(1, size(a, 2));
  b(col) = row;
  b = b + 1;
end