function k_nn_pixels()
  % k-nearest neighbours in the normalized image coordinates
  % 
  % Pavel Trutman
  % INRIA, 2016

  deepagpaths;

  % prepare data
  load('../../data/paris/correspondences.mat');
  
  % properties
  margin = 1/7;
  batchSize = 10;
  
  % filter cams that their correspondences are close to each other
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
  
  %load('pixels_syntetic.mat');
  
  fprintf(['Size of training dataset: ', num2str(size(u_tr, 2)), '.\n']);
  fprintf(['Size of validating dataset: ', num2str(size(u_val, 2)), '.\n']);
  
  % nearest neighbour for points validating dataset
  fprintf('Searching for nearest neighbour for points in validating dataset.\n'); tic;
  
  u1x_tr = u_tr(1:7, :);
  u1y_tr = u_tr(8:14, :);
  u2x_tr = u_tr(15:21, :);
  u2y_tr = u_tr(22:28, :);
  
  % initialize the parallel pool and the wrappers
  pool = gcp();
  u1x_trWrap = WorkerObjWrapper(u1x_tr);
  u1y_trWrap = WorkerObjWrapper(u1y_tr);
  u2x_trWrap = WorkerObjWrapper(u2x_tr);
  u2y_trWrap = WorkerObjWrapper(u2y_tr);

  iterations = idivide(int32(size(u_val, 2)), batchSize, 'ceil');
  results(iterations) = parallel.FevalFuture();
  
  % submit data to parallel computation in batches
  repS = {};
  for i = 1:iterations
    iLow = (i - 1)*batchSize + 1;
    iTop = min(size(u_val, 2), i*batchSize);
    repS = adprintf(repS, [num2str(iLow), '/', num2str(size(u_val, 2))]);
    u1x_val = u_val(1:7, iLow:iTop);
    u1y_val = u_val(8:14, iLow:iTop);
    u2x_val = u_val(15:21, iLow:iTop);
    u2y_val = u_val(22:28, iLow:iTop);
    
    results(i) = parfeval(pool, @computeDist, 2, u1x_val, u1y_val, u2x_val, u2y_val, u1x_trWrap, u1y_trWrap, u2x_trWrap, u2y_trWrap, margin);
    repS = rmprintf(repS);
  end
  
  % retrieve results from parallel computation
  distances = NaN(size(u_val, 2), 1);
  indices = NaN(size(u_val, 2), 1);
  repS = {};
  for i = 1:iterations
    t = tic;
    [ii, dist, ind] = fetchNext(results);
    if ii ~= iterations
      distances(((ii-1)*batchSize+1):(ii*batchSize)) = dist;
      indices(((ii-1)*batchSize+1):(ii*batchSize)) = ind;
    else
      distances(((ii-1)*batchSize+1):end) = dist;
      indices(((ii-1)*batchSize+1):end) = ind;
    end
    repS = rmprintf(repS);
    repS = adprintf(repS, sprintf([num2str((i - 1)*batchSize + 1), '/', num2str(size(u_val, 2)), ': %1.3f s'], toc(t)));
  end
  rmprintf(repS);
  
  % clean up
  fprintf(' ... done %s.\n', sec2hms(toc));
  delete(u1x_trWrap);
  delete(u1y_trWrap);
  delete(u2x_trWrap);
  delete(u2y_trWrap);
end

%% Functions
  
function cix = findGoodCams(u, margin)
  % Finds cameras from dataset with correspondences that are far away from
  % each other.
  % 
  % u = coordinates of correspondences, size: 28xdatasetSize, each column
  %     is one pair of camerase in form: [u1x; u1y; u2x; u2y]
  % margin = distance how near two points in one image can be
  % cix = logical array of dataset size selecting good cameras
  
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
    
    % check constraint for the first cam
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
    
    % check constraint for the second cam
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

function [dist, index] = computeDist(u1x_val, u1y_val, u2x_val, u2y_val, u1x_trWrap, u1y_trWrap, u2x_trWrap, u2y_trWrap, margin)
  % Computes minimal distance for each point in validating dataset from its
  % nearest neighbour. Resulting distance is maximal euclidean distance of
  % each correspondence.
  % 
  % This function is used for parallel computation, therefore subset of
  % validating dataset of size batchSize is selected to minimize overhead
  % of function call.
  % 
  % u1x_val = ...
  % u1y_val = ...
  % u2x_val = ...
  % u2y_val = coordinates of the correspondences, size 7xbatchSize
  % u1x_trWrap = ...
  % u1y_trWrap = ...
  % u2x_trWrap = ...
  % u2y_trWrap = wrapper containing coordinates of training dataset
  %              correspondences, size: 7xdatasetSize, wrapper is used to
  %              minimize amount of data sent to parallel workers
  % margin = maximal distance of coordinate of validating point from
  %          training point
  % dist = distance of validating point from training point, NaN if margin
  %        constraint violated, size: batchSize
  % index = index of nearest neighbour in training dataset, NaN if margin
  %         constraint violated, size: batchSize
  
  dist = NaN(size(u1x_val, 2), 1);
  index = NaN(size(u1x_val, 2), 1);

  % unwrap the wrappers
  u1x_tr = u1x_trWrap.Value;
  u1y_tr = u1y_trWrap.Value;
  u2x_tr = u2x_trWrap.Value;
  u2y_tr = u2y_trWrap.Value;
  
  % for each validating point in batch
  for i = 1:size(u1x_val, 2)
    % check the margin constraint, not swapped 1-1
    ordernoswap = zeros(7, size(u1x_tr, 2));
    noswap = true(1, size(u1x_tr, 2));
    for j = 1:7
      hxBound = repmat(u1x_val(j, i) + margin, 7, sum(noswap));
      hyBound = repmat(u1y_val(j, i) + margin, 7, sum(noswap));
      lxBound = repmat(u1x_val(j, i) - margin, 7, sum(noswap));
      lyBound = repmat(u1y_val(j, i) - margin, 7, sum(noswap));
      match = ((u1x_tr(:, noswap) < hxBound) & (u1x_tr(:, noswap) > lxBound) & (u1y_tr(:, noswap) < hyBound) & (u1y_tr(:, noswap) > lyBound));
      ordernoswap(j, noswap) = findFirstInCol(match);
      noswap(noswap) = (sum(match, 1) == 1);
    end
    % check uniqueness of the order
    noswap(noswap) = all(sort(ordernoswap(:, noswap)) == repmat((1:7)', 1, sum(noswap)));
    % check the margin constraint, not swapped 2-2
    idx = (ordernoswap(:, noswap) - 1)'*sum(noswap) + ndgrid(1:sum(noswap), 1:7);
    u2x_sel = u2x_tr(:, noswap)';
    u2x_sel(:) = u2x_sel(idx);
    u2y_sel = u2y_tr(:, noswap)';
    u2y_sel(:) = u2y_sel(idx);
    hxBound = repmat(u2x_val(:, i) + margin, 1, sum(noswap));
    hyBound = repmat(u2y_val(:, i) + margin, 1, sum(noswap));
    lxBound = repmat(u2x_val(:, i) - margin, 1, sum(noswap));
    lyBound = repmat(u2y_val(:, i) - margin, 1, sum(noswap));
    match = all((u2x_sel' < hxBound) & (u2x_sel' > lxBound) & (u2y_sel' < hyBound) & (u2y_sel' > lyBound));
    noswap(noswap) = match;
    
    % check the margin constraint, swapped 1-2
    orderswap = zeros(7, size(u1x_tr, 2));
    swap = true(1, size(u1x_tr, 2));
    for j = 1:7
      hxBound = repmat(u1x_val(j, i) + margin, 7, sum(swap));
      hyBound = repmat(u1y_val(j, i) + margin, 7, sum(swap));
      lxBound = repmat(u1x_val(j, i) - margin, 7, sum(swap));
      lyBound = repmat(u1y_val(j, i) - margin, 7, sum(swap));
      match = ((u2x_tr(:, swap) < hxBound) & (u2x_tr(:, swap) > lxBound) & (u2y_tr(:, swap) < hyBound) & (u2y_tr(:, swap) > lyBound));
      orderswap(j, swap) = findFirstInCol(match);
      swap(swap) = (sum(match, 1) == 1);
    end
    % check uniqueness of the order
    swap(swap) = all(sort(orderswap(:, swap)) == repmat((1:7)', 1, sum(swap)));
    % check the margin constraint, swapped 2-1
    idx = (orderswap(:, swap) - 1)'*sum(swap) + ndgrid(1:sum(swap), 1:7);
    u1x_sel = u1x_tr(:, swap)';
    u1x_sel(:) = u1x_sel(idx);
    u1y_sel = u1y_tr(:, swap)';
    u1y_sel(:) = u1y_sel(idx);
    hxBound = repmat(u2x_val(:, i) + margin, 1, sum(swap));
    hyBound = repmat(u2y_val(:, i) + margin, 1, sum(swap));
    lxBound = repmat(u2x_val(:, i) - margin, 1, sum(swap));
    lyBound = repmat(u2y_val(:, i) - margin, 1, sum(swap));
    match = all((u1x_sel' < hxBound) & (u1x_sel' > lxBound) & (u1y_sel' < hyBound) & (u1y_sel' > lyBound));
    swap(swap) = match;
    
    % compute the distances for each training point, that does not violate
    % the margin constraint
    noswapIdx = find(noswap);
    noswapDist = zeros(size(noswapIdx, 2), 1);
    for j = noswapIdx
      idx = ordernoswap(:, j);
      noswapDist(j) = max([sum([(u1x_val(:, i) - u1x_tr(idx, j))'; (u1y_val(:, i) - u1y_tr(idx, j))'].^2) sum([(u2x_val(:, i) - u2x_tr(idx, j))'; (u2y_val(:, i) - u2y_tr(idx, j))'].^2)]);
    end
    swapIdx = find(swap);
    swapDist = zeros(size(swapIdx, 2), 1);
    for j = swapIdx
      idx = orderswap(:, j);
      swapDist(j) = max([sum([(u1x_val(:, i) - u2x_tr(idx, j))'; (u1y_val(:, i) - u2y_tr(idx, j))'].^2) sum([(u2x_val(:, i) - u1x_tr(idx, j))'; (u2y_val(:, i) - u1y_tr(idx, j))'].^2)]);
    end
    
    % select the minimal distance and remember the index
    allIdx = [noswapIdx, swapIdx];
    if size(allIdx, 2) > 0
      allDist = [noswapDist; swapDist];
      [minDist, idx] = min(allDist);
      dist(i) = minDist;
      index(i) = allIdx(idx);
    end
    
  end
  
  function b = findFirstInCol(a)
    % finds first nonzero element in each column
    % 
    % a = matrix to search in, size: rowsxcols
    % b = vector of indices of first nonzero element in each column, size:
    %     1xcols
    
    a = cumsum(a);
    [row, col] = find(a == 0);
    b = zeros(1, size(a, 2));
    b(col) = row;
    b = b + 1;
  end

end