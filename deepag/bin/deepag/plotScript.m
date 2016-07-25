% Plot script
% 
% Pavel Trutman
% INRIA, 2016

file = 'nn_4L_GPU-160720-152910.mat';
load(['results/' file]);

if ~exist('plotPerEpoch', 'var')
  plotPerEpoch = 1;
end

figure;
hold on;

x = (0:size(tr_error, 2)-1)*plotPerEpoch;

semilogy(x, tr_error);
semilogy(x, val_error);
legend('Training dataset', 'Validating dataset', 'Location', 'southoutside');
xlabel('Epoch');
ylabel('Error');
xlim([0 x(end)]);
ylim([min([tr_error val_error]) max([tr_error val_error])]);
grid on;
hold off;