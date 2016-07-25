% Plot k-nn script
% 
% Pavel Trutman
% INRIA, 2016

file = 'k_nn_GPU-160725-150925.mat';
load(['results/' file]);

figure;
histogram(val_error, 'Normalization', 'probability');
xlabel('Error [px]');
ylabel('Probability');
title('Histogram of errors.');

if exist('val_error_k', 'var')
  figure;
  hold on;
  clear legendInfo;
  for i = 1:k
    histogram(val_error_k(i, :), 'Normalization', 'probability');
    legendInfo{i} = (['Error from ', num2str(i), '. neighbour']);
  end
  xlabel('Error [px]');
  ylabel('Probability');
  legend(legendInfo);
  title('Histogram of errors to each neighbour.');
  hold off;

  figure;
  histogram(std(val_error_k), 'Normalization', 'probability');
  xlabel('Standard deviation of error [px]');
  ylabel('Probability');
  title('Histogram of standard deviation of error.');
end
