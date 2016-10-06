
val_file='../../data/paris/features_F_synth.mat';
tr_file='../../data/paris/features_F_synth_sample_1K.mat';
run('test')

val_file='../../data/paris/features_F_nsynthrep.mat';
tr_file='../../data/paris/features_F_synth_sample_1K.mat';
run('test')

val_file='../../data/paris/features_F_synth.mat';
tr_file='../../data/paris/features_F_nsynthrep_sample_1K.mat';
run('test')

val_file='../../data/paris/features_F_nsynthrep.mat';
tr_file='../../data/paris/features_F_nsynthrep_sample_1K.mat';
run('test')
