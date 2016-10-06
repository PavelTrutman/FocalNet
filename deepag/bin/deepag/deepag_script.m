% Automatic script for data preparation
% 
% Pavel Trutman
% INRIA, 2016

deepag_parallel();
ps.Load = 1;
ps.Normalize = 1;
ps.RemoveOutliers = 1;
ps.Divide = 1;
ps.Prepare = 1;
ps.quietMode = 1;

deepag_parallel(ps);
exit;