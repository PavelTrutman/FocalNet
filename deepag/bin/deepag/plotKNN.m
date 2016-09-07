% Plot k-nn script
% 
% Pavel Trutman
% INRIA, 2016

% load file
file = 'k_nn-160906-142446-1K-noise in the u.mat';
load(['results/' file]);

% plot histogram of error
figure;
histogram(val_error, 'Normalization', 'probability');
xlabel('Error [px]');
ylabel('Probability');
title('Histogram of errors.');

% plot histogram of error from each neighbour
if exist('val_error_k', 'var')
  figure;
  hold on;
  clear legendInfo;
  for i = 1:k
    histogram(val_error_k(i, :), 'Normalization', 'probability');
    legendInfo{i} = (['Error from ', num2str(i), '. neighbour']); %#ok<SAGROW>
  end
  xlabel('Error [px]');
  ylabel('Probability');
  legend(legendInfo);
  title('Histogram of errors to each neighbour.');
  hold off;
end
mean(val_error)