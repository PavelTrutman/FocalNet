% deepag.m - Deep Learning 4 Algebraic Geometry
% T. Pajdla, pajdla@gmail.cz
% 1 Jan 2016 - 31 Dec 2016
%
function deepag(pst, pdAllt, datasett)
  persistent ps;
  persistent DEPPAGRunning;
  persistent pdAll;
  persistent dataset;
  if exist('pst', 'var')
    ps = pst;
  end
  if exist('pdAllt', 'var')
    pdAll = pdAllt;
  end
  if exist('datasett', 'var')
    dataset = datasett;
  end

  %% GUI structures
  %% Clear, close all if for the first time from the command line.
  if isempty(DEPPAGRunning) || (DEPPAGRunning ~= true)
    deepagpaths;          % set paths to libraries
    close all             % hardcore close all
    %clear                 % clear variables
    %% Initialize parameters structure ps
    ps.Data           = 'paris'; % Data path
    ps.Load           = false;   % Load data
    ps.ShowMap        = false;   % show map
    ps.Normalize      = false;   % Normalize
    ps.RemoveOutliers = false;   % Remove focal length outliers
    ps.Divide         = false;   % Divide points into train/val/test subsets
    ps.Hist           = false;   % Show histograms of the datasets
    ps.Prepare        = false;   % Choose a pair of cameras with many common points  at random and choose a set of 7 points at random
    ps.Exit           = false;   % exit & close all
    %% Costruct GUI for the ps parameters
    gf(1) = subfig(4,9,1); % GUI window
    % gf(1) = subfig([1 50 763 127],figure); % GUI window
    set(gcf,'closerequestfcn',[]);     % make GUI unclosable
    ps = pargui(ps, [size(fieldnames(ps),1)+1 1], 'deepag');   % create the GUI 
    parcutscreen(gf(1));               % shrink the screen
    clear gf
    %
    DEPPAGRunning = true;                  % running    
  else %close all except the unclosable GUI and keep existing parameters
    deepagInit
  end
  
  %% Execs ------------------------------------------------------------------
  %% Close and clear all and exit
  if ps.Exit, deepagExit; return; end
  
  %% Load data
  if ps.Load
    fprintf('%s: Load data: ',ps.Data); tic;
    repS = adprintf({},'Loading '); % progress report
    switch pc.dataType
      case 'paris'
        for i = 1:4
          fn = fullfile(pc.dataPath, ps.Data,sprintf('RS_Paris_%d.mat', i));
          repS = adprintf(repS, sprintf('%s %s', sec2hms(toc), fn));
          d = load(fn);
          d = d.RS_part;
          pd.X{i} = d.X.pts;
          pd.C{i}.cix = d.cams;
          pd.C{i}.cam = d.cam;
          repS = rmprintf(repS);
        end
        is = load(fullfile(pc.dataPath,ps.Data,'paris_imres.mat'));
        pd.is = is.imres;
                
        % merge together
        pdAll.X = [pd.X{1}, pd.X{2}, pd.X{3}, pd.X{4}];
        pdAll.C.cix = [pd.C{1}.cix, pd.C{2}.cix, pd.C{3}.cix, pd.C{4}.cix];
        pdAll.C.cam = [pd.C{1}.cam, pd.C{2}.cam, pd.C{3}.cam, pd.C{4}.cam];
        % one cam is missing
        pdAll.C.cam = {pdAll.C.cam{1:16892-1} pdAll.C.cam{570} pdAll.C.cam{16892:end}};
        pdAll.is = pd.is;
        clear pd;
                
      otherwise
        repS = adprintf(repS,'Unknown pc.dataType');
    end
    repS = rmprintf(repS);
    repS = adprintf(repS,sprintf(' Loaded %s.',sec2hms(toc)));
    rmprintf(repS);
    clear fn repS d i is d;
    fprintf(' ... done %s.\n',sec2hms(toc));
    assignin('base', 'pdAll', pdAll);
  end
  
  %% Show map of 3D points and cameras' centers
  if ps.ShowMap
    fprintf('%s: Showing Map: ',ps.Data); tic;
    subfig(2,2,1); hold on
    plot3d(pdAll.X, ['.' pc.clrs(1)]);
    axis equal; grid;
    title('3D points');
    subfig(2,2,2); hold on
    T = cellfunuo(@(c)null(c(1).P), pdAll.C.cam(pdAll.C.cix));
    T = cat(2,T{:}); T = h2a(T);
    plot3d(T,['.' pc.clrs(1)]);
    axis equal; grid
    title('Camera centers')
    subfig(2,2,3); hold on
    plot3d(pdAll.X, ['.' pc.clrs(1)]);
    plot3d(T,['.' pc.clrs(2)]);
    axis equal; grid;
    title('3D points and camera centers');
    fprintf(' ... done %s.\n',sec2hms(toc));
    clear T
  end
  
  %% Normalize data from pixels to image widths
  if ps.Normalize
    fprintf('%s: Normalizing data: ', ps.Data); tic;
    for i = pdAll.C.cix
      pdAll.C.cam{i}.normFactor = max(pdAll.is(i, :))/2;
      pdAll.C.cam{i}.u = [pdAll.C.cam{i}.u(1, :)./pdAll.C.cam{i}.normFactor; pdAll.C.cam{i}.u(2, :)./pdAll.C.cam{i}.normFactor];
      pdAll.C.cam{i}.f = pdAll.C.cam{i}.f/pdAll.C.cam{i}.normFactor;
    end
    fprintf(' ... done %s.\n',sec2hms(toc));
    clear i;
    assignin('base', 'pdAll', pdAll);
  end
  
  %% Remove cameras with extra large focal length
  if ps.RemoveOutliers
    fprintf('%s: Removing outliers: ', ps.Data); tic;
    fprintf('\n');
    
    f = zeros(1, size(pdAll.C.cam, 2));

    for i = pdAll.C.cix
      f(i) = pdAll.C.cam{i}.f;
    end
    inliers = pdAll.C.cix;
    
    inliersLast = size(inliers, 2) + 1;
    while (inliersLast - size(inliers, 2)) ~= 0
      inliersLast = size(inliers, 2);
      outliers = find(f > (mean(f(inliers)) + 4*std(f(inliers))));
      inliers = setdiff(inliers, outliers);
    end
    pdAll.C.inliers = inliers;
    fprintf(['Removed ', num2str(size(outliers, 2)), 'cams.\n']);
    fprintf(' ... done %s.\n',sec2hms(toc));
    clear f i outliers inliers inliersLast;
    assignin('base', 'pdAll', pdAll);
  end
  
  %% Divide 3D points and cams into three non-intersecting datasets
  if ps.Divide
    fprintf('%s: Showing Divided Map: \n', ps.Data); tic;
    fn = fullfile(pc.dataPath,ps.Data, 'boundary.mat');
    load(fn);
    subfig(1,1,1); hold on
    x = pdAll.X([1 3], :);
    plot(x(1, :), x(2, :), ['.' pc.clrs(1)]);
    axis equal; grid;
    title('Division into training, validating and testing sets.');
    hold on;
    clear fn;
    
    % divide 3D points
    dataset = cell(1, 3);
    for i  = 1:3
      plot(boundary{i}(1, :), boundary{i}(2, :), ['-' pc.clrs(i + 1)]); %#ok<USENS>
      dataset{i}.Xmask = find(inpolygon(x(1, :), x(2, :), boundary{i}(1, :), boundary{i}(2, :)));
      fprintf(['Points in dataset', int2str(i), ': ', int2str(size(dataset{i}.Xmask, 2)), '\n']);
    end
    drawnow;
    clear i x boundary;
    
    % divide cameras
    camDatasetIdx = zeros(1, size(pdAll.C.cam, 2));
    camDataset = repmat(struct('Xmask', [], 'umask', [], 'cix', 0), 1, size(pdAll.C.cam, 2));
    for i = pdAll.C.inliers
      repS = adprintf({}, ['Parsing camera ' int2str(i), '/', int2str(size(pdAll.C.cam, 2))]);
      camXmask = pdAll.C.cam{i}.Xmask;
      
      Xmask = cell(1, 3);
      idx = cell(1, 3);
      [Xmask{1}, idx{1}, ~] = intersect(camXmask, dataset{1}.Xmask);
      [Xmask{2}, idx{2}, ~] = intersect(camXmask, dataset{2}.Xmask);
      [Xmask{3}, idx{3}, ~] = intersect(camXmask, dataset{3}.Xmask);
      [~, idataset] = max([size(Xmask{1}, 2), size(Xmask{2}, 2), size(Xmask{3}, 2)]);
      
      % skip cameras that can see points from other datasets
      if (1 ~= idataset) && (size(Xmask{1}, 2) > 0)
        rmprintf(repS);
        continue;
      end
      if (2 ~= idataset) && (size(Xmask{2}, 2) > 0)
        rmprintf(repS);
        continue;
      end
      if (3 ~= idataset) && (size(Xmask{3}, 2) > 0)
        rmprintf(repS);
        continue;
      end
      
      camDataset(i).Xmask = Xmask{idataset};
      camDataset(i).umask = idx{idataset}';
      camDataset(i).cix = i;
      camDatasetIdx(i) = idataset;
      rmprintf(repS);
    end
    
    clear inliers cams;
    for i = 1:3
      dataset{i}.cam = num2cell(camDataset(camDatasetIdx == i));
      fprintf(['Cameras in dataset', int2str(i), ': ', int2str(size(dataset{i}.cam, 2)), '\n']);
    end
    fprintf(' ... done %s.\n',sec2hms(toc));
    clear repS i j cam Xmask idx idataset camDataset camDatasetIdx;
    assignin('base', 'dataset', dataset);
  end
  
  %% Plot many histograms
  if ps.Hist
    fprintf('%s: Creating histograms: ', ps.Data); tic;
    h = cell(1, 3);
    for i = 1:3
      h{i}.f = zeros(1, size(dataset{i}.cam, 2));
      h{i}.fNorm = zeros(1, size(dataset{i}.cam, 2));
      h{i}.width = zeros(1, size(dataset{i}.cam, 2));
      h{i}.height = zeros(1, size(dataset{i}.cam, 2));
      h{i}.ratio = zeros(1, size(dataset{i}.cam, 2));
      
      % prepare data
      for j = 1:size(dataset{i}.cam, 2)
        h{i}.f(j) = pdAll.C.cam{dataset{i}.cam{j}.cix}.f*pdAll.C.cam{dataset{i}.cam{j}.cix}.normFactor;
        h{i}.fNorm(j) = pdAll.C.cam{dataset{i}.cam{j}.cix}.f;
        h{i}.width(j) = max(pdAll.is(dataset{i}.cam{j}.cix, :));
        h{i}.height(j) = min(pdAll.is(dataset{i}.cam{j}.cix, :));
        h{i}.ratio(j) = h{i}.width(j)/h{i}.height(j);
      end
    end

    % focal length histogram
    figure;
    hold on;
    for i = 1:3
      histogram(h{i}.f, 'Normalization', 'probability', 'BinWidth', 20, 'FaceAlpha', 0.5);
      xlabel('Focal length [px]');
      ylabel('Probability');
    end
    legend('Testing dataset', 'Validating dataset', 'Training dataset');
    %xlim([0 10000]);
    hold off;
    
    % normalized focal length histogram
    figure;
    hold on;
    for i = 1:3
      histogram(h{i}.fNorm, 'Normalization', 'probability', 'BinWidth', 0.1, 'FaceAlpha', 0.5);
      xlabel('Normalized focal length [width of image]');
      ylabel('Probability');
    end
    legend('Testing dataset', 'Validating dataset', 'Training dataset');
    xlim([0 10]);
    hold off;

    % width histogram
    figure;
    hold on;
    for i = 1:3
      histogram(h{i}.width, 'Normalization', 'probability', 'BinWidth', 10, 'FaceAlpha', 0.5);
      xlabel('Image width [px]');
      ylabel('Probability');
    end
    legend('Testing dataset', 'Validating dataset', 'Training dataset');
    hold off;

    % height histogram
    figure;
    hold on;
    for i = 1:3
      histogram(h{i}.height, 'Normalization', 'probability', 'BinWidth', 10, 'FaceAlpha', 0.5);
      xlabel('Image height [px]');
      ylabel('Probability');
    end
    legend('Testing dataset', 'Validating dataset', 'Training dataset');
    hold off;

    % ratio histogram
    figure;
    hold on;
    for i = 1:3
      histogram(h{i}.ratio, 'Normalization', 'probability', 'BinWidth', 0.01, 'FaceAlpha', 0.5);
      xlabel('Image ratio (width/height)');
      ylabel('Probability');
    end
    legend('Testing dataset', 'Validating dataset', 'Training dataset');
    xlim([1 2]);
    hold off;
    
    assignin('base', 'h', h);
    fprintf(' ... done %s.\n',sec2hms(toc));
    clear i j;
  end
  
  %% Randomly sample the datasets and generate feature vectors and correspondences
  if ps.Prepare
    cameraPairsNum = 10000;
    minPointsInCommon = 15;
    pointsNum = 7;
    perCameraPair = 100;
    coefsSize = 245;
    fSize = 2;
    normSize = 2;
    features = cell(1, 3);
    correspondences = cell(1, 3);
    
    fprintf('%s: Preparing data: \n',ps.Data); tic;
    cams = cellfun(@camU, pdAll.C.cam, 'UniformOutput', false);
    
    % for each dataset
    for i = 1:3
      repS = adprintf({}, ['Preparing dataset ' int2str(i), ': ']);
      features{i}.coefs = zeros(coefsSize, cameraPairsNum*perCameraPair);
      features{i}.f = zeros(fSize, cameraPairsNum*perCameraPair);
      features{i}.norm = zeros(normSize, cameraPairsNum*perCameraPair);
      correspondences{i}.u = zeros(pointsNum*4, cameraPairsNum*perCameraPair);
      correspondences{i}.f = zeros(fSize, cameraPairsNum*perCameraPair);
      correspondences{i}.norm = zeros(normSize, cameraPairsNum*perCameraPair);
      dataseti = dataset{i}.cam;
      
      % randomly sample camera pairs
      j = 1;
      while j <= cameraPairsNum
        repS = adprintf(repS, ['camera pair ', int2str(j), '/', int2str(cameraPairsNum)]);
        r = randperm(size(dataset{i}.cam, 2), 2);
  
        [coefs, u] = prepare(r, dataseti, cams, perCameraPair, pointsNum, minPointsInCommon, coefsSize);
        if any(isnan(coefs))
          repS = rmprintf(repS);
          continue;
        end
        
        features{i}.coefs(:, ((j-1)*perCameraPair+1):(j*perCameraPair)) = coefs;
        features{i}.f(:, ((j-1)*perCameraPair+1):(j*perCameraPair)) = repmat([pdAll.C.cam{dataset{i}.cam{r(1)}.cix}.f; pdAll.C.cam{dataset{i}.cam{r(2)}.cix}.f], 1, perCameraPair);
        features{i}.norm(:, ((j-1)*perCameraPair+1):(j*perCameraPair)) = repmat([pdAll.C.cam{dataset{i}.cam{r(1)}.cix}.normFactor; pdAll.C.cam{dataset{i}.cam{r(2)}.cix}.normFactor], 1, perCameraPair);
        correspondences{i}.u(:, ((j-1)*perCameraPair+1):(j*perCameraPair)) = u;
        correspondences{i}.f(:, ((j-1)*perCameraPair+1):(j*perCameraPair)) = repmat([pdAll.C.cam{dataset{i}.cam{r(1)}.cix}.f; pdAll.C.cam{dataset{i}.cam{r(2)}.cix}.f], 1, perCameraPair);
        correspondences{i}.norm(:, ((j-1)*perCameraPair+1):(j*perCameraPair)) = repmat([pdAll.C.cam{dataset{i}.cam{r(1)}.cix}.normFactor; pdAll.C.cam{dataset{i}.cam{r(2)}.cix}.normFactor], 1, perCameraPair);
        
        repS = rmprintf(repS);
        j = j + 1;
      end
      rmprintf(repS);
      
    end
    
    % remove constant elements
    X_std = std([features{1}.coefs, features{2}.coefs features{3}.coefs], [], 2);
    for i = 1:3
      features{i}.coefs = features{i}.coefs(X_std ~= 0, :);
    end
    
    clear repS i j k l cameraPairsNum minPointsInCommon pointsNum perCameraPair coefsSize fSize normSize r XmaskCommon idx1 idx2 rk coefs X_std keep;
    
    assignin('base', 'features', features);
    assignin('base', 'correspondences', correspondences);
    
    % save feature vectors into file
    fn = fullfile(pc.dataPath,ps.Data, 'features.mat');
    fprintf(['Saving feature vectors into ', fn, '\n']);
    tr_coefs = features{3}.coefs; %#ok<NASGU>
    save(fn, 'tr_coefs', '-v7.3');
    clear tr_coefs;
    val_coefs = features{2}.coefs; %#ok<NASGU>
    save(fn, 'val_coefs', '-append');
    clear val_coefs;
    tst_coefs = features{1}.coefs; %#ok<NASGU>
    save(fn, 'tst_coefs', '-append');
    clear tst_coefs;
    tr_f = features{3}.f; %#ok<NASGU>
    save(fn, 'tr_f', '-append');
    clear tr_f;
    val_f = features{2}.f; %#ok<NASGU>
    save(fn, 'val_f', '-append');
    clear val_f;
    tst_f = features{1}.f; %#ok<NASGU>
    save(fn, 'tst_f', '-append');
    clear tst_f;
    tr_norm = features{3}.norm; %#ok<NASGU>
    save(fn, 'tr_norm', '-append');
    clear tr_norm;
    val_norm = features{2}.norm; %#ok<NASGU>
    save(fn, 'val_norm', '-append');
    clear val_norm;
    tst_norm = features{1}.norm; %#ok<NASGU>
    save(fn, 'tst_norm', '-append');
    clear tst_norm;
    
    % save correspondences into file
    fn = fullfile(pc.dataPath,ps.Data, 'correspondences.mat');
    fprintf(['Saving correspondences into ', fn, '\n']);
    corr_tr = correspondences{3}; %#ok<NASGU>
    corr_val = correspondences{2}; %#ok<NASGU>
    corr_tst = correspondences{1}; %#ok<NASGU>
    save(fn, 'corr_tr', 'corr_val', 'corr_tst',  '-v7.3');
    clear corr_tr corr_val corr_tst;
    
    fprintf(' ... done %s.\n',sec2hms(toc));
    
  end

  assignin('base', 'ps', ps);
  assignin('base', 'pdAll', pdAll);
  assignin('base', 'dataset', dataset);
  
end

%% Functions

function [u] = camU(cam)
  % returns coordinates of correspondences from cam
  % 
  % cam = structure representing camera
  % u = coordinates of correspondences

  if isempty(cam)
    u = [];
  else
    u = cam.u;
  end
  
end

function [coefsOut, uOut] = prepare(r, dataset, cams, perCameraPair, pointsNum, minPointsInCommon, coefsSize)
  % For two cameras returns feature vector and coordinates of
  % correspondences.
  % 
  % r = [index of the first camera, index of the second camera]
  % dataset = cameras from dataset
  % cams = all cameras
  % perCameraPair = number of samples for each camera pair
  % pointsNum = how many correspondences per sample
  % minPointsInCommon = minimal number of points that the two cameras have
  %                     to have in common, is used to prevent to sample the
  %                     same correspondences
  % coefsSize = expected size feature vector
  % coefsOut = features vectors, size: coefsSize x perCameraPair
  % uOut = coordinates of correspondences, size: 4*pointsNum x perCameraPair
  
  % prepare data
  datasetCam1 = dataset{r(1)};
  datasetCam2 = dataset{r(2)};
  u1 = cams{datasetCam1.cix};
  u2 = cams{datasetCam2.cix};

  % compute 3D points that can see both cameras
  [XmaskCommon, idx1, idx2] = intersect(datasetCam1.Xmask, datasetCam2.Xmask);
  if size(XmaskCommon, 2) < minPointsInCommon
    coefsOut = NaN;
    uOut = NaN;
    return;
  end

  % randomly sample correspondences
  k = 1;
  coefsOut = zeros(coefsSize, perCameraPair);
  uOut = zeros(pointsNum*4, perCameraPair);
  while k <= perCameraPair
    rk = randperm(size(XmaskCommon, 2), pointsNum);
    coefs = get_features(u1(:, datasetCam1.umask(idx1(rk)))', u2(:, datasetCam2.umask(idx2(rk)))');
    if any(isnan(coefs))
      continue;
    end
    coefsOut(:, k) = coefs;
    uOut(:, k) = [reshape(u1(:, datasetCam1.umask(idx1(rk))), [], 1); reshape(u2(:, datasetCam2.umask(idx2(rk))), [], 1)];
    k = k + 1;
  end
end