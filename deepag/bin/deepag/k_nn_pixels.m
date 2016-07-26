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
  for i = 1:size(u_val, 2)
    u1_val = [u_val(1:7, i)'; u_val(8:14, i)'];
    u2_val = [u_val(15:21, i)'; u_val(22:28, i)'];
    
    dist = Inf(size(u_tr, 2), 1);
    for k = 1:size(u_tr, 2)
      u1_tr = [u_tr(1:7, i)'; u_tr(8:14, i)'];
      u2_tr = [u_tr(15:21, i)'; u_tr(22:28, i)'];
      
      order1 = zeros(7, 1);
      swap = false;
      skip = false;
      for j = 1:7
        hBound = repmat(u1_val(:, j) + repmat(margin, 2, 1), 1, 7);
        lBound = repmat(u1_val(:, j) - repmat(margin, 2, 1), 1, 7);
        fix = find(all((u1_tr < hBound) & (u1_tr > lBound)));
        
        if size(fix, 2) ~= 1
          swap = true;
          break;
        end
        order1(j) = fix;
      end
      
      if ~swap
        order2 = zeros(7, 1);
        for j = 1:7
          hBound = repmat(u2_val(:, j) + repmat(margin, 2, 1), 1, 7);
          lBound = repmat(u2_val(:, j) - repmat(margin, 2, 1), 1, 7);
          fix = find(all((u2_tr < hBound) & (u2_tr > lBound)));
        
          if size(fix, 2) ~= 1
            swap = true;
            break;
          end
          order2(j) = fix;
        end
      end
      
      if swap
        % swap the cams
        
        order1 = zeros(7, 1);
        for j = 1:7
          hBound = repmat(u1_val(:, j) + repmat(margin, 2, 1), 1, 7);
          lBound = repmat(u1_val(:, j) - repmat(margin, 2, 1), 1, 7);
          fix = find(all((u2_tr < hBound) & (u2_tr > lBound)));
        
          if size(fix, 2) ~= 1
            skip = true;
            break;
          end
          order1(j) = fix;
        end
        
        if skip
          continue;
        end
        
        order2 = zeros(7, 1);
        for j = 1:7
          hBound = repmat(u2_val(:, j) + repmat(margin, 2, 1), 1, 7);
          lBound = repmat(u2_val(:, j) - repmat(margin, 2, 1), 1, 7);
          fix = find(all((u1_tr < hBound) & (u1_tr > lBound)));
        
          if size(fix, 2) ~= 1
            skip = true;
            break;
          end
          order2(j) = fix;
        end
        
        if skip
          break;
        end
        
      end
      
      % points are close
      if swap
        [u2_tr, u1_tr] = deal(u1_tr, u2_tr);
      end
        
      dist(k) = max([sum((u1_tr(:, order1) - u1_val).^2, 2), sum((u2_tr(:, order2) - u2_val).^2, 2)]);
      
    end
    [m, idx] = min(dist);
    fprintf([num2str(i), '/', num2str(size(u_val, 2)), ': ', num2str(m), '\n']);
    
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