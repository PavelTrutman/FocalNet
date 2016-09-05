% Script for generating synthetic Fundamental matrix feature vectors
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016
function F_features_u()

trSize = 1*1000;
valSize = 1*1000;

tr_coefs=zeros(9,trSize);
tr_f=zeros(2,trSize);
tr_norm=zeros(2,trSize);
for i=1:trSize
    Fg=getRandomFg();
    fg = F2f1f2(Fg);
    tr_coefs(:,i)=reshape(Fg,9,1);
    tr_f(:,i)=fg;
    tr_norm(:,i)=[1 1];
end
save('../../data/paris/features_F_synth.mat', 'tr_coefs', 'tr_norm', 'tr_f', '-v7.3');
clear tr_f tr_coefs tr_norm;

val_coefs=zeros(9,valSize);
val_f=zeros(2,valSize);
val_norm=zeros(2,valSize);
for i=1:valSize
    Fg=getRandomFg();
    fg = F2f1f2(Fg);
    val_coefs(:,i)=reshape(Fg,9,1);
    val_f(:,i)=fg;
    val_norm(:,i)=[1 1];
end
save('../../data/paris/features_F_synth_sample_10K.mat', 'val_coefs', 'val_norm', 'val_f', '-v7.3');
clear val_f val_coefs val_norm;
end

function Fg=getRandomFg()
minlen=0; diflen=1;
noise=true;
% focal lengths
f1 = minlen+diflen*rand(1,1);
f2 = minlen+diflen*rand(1,1);
sc = max(f1,f2);
% Internal camera calibration
K1 = [f1  0    0
    0   f1   0
    0   0    1];
K2 = [f2 0  0
    0  f2 0
    0  0  1];
% Rotations
R1 = eye(3);
R2 = a2r(rand(3,1),pi/10);
% Projection matrices
P1 = K1*R1*[eye(3)     [-sc;0;0]];
P2 = K2*R2*[eye(3)  sc*[rand(2,1);0]];
% Fundamental matrix
Fg = PP2F(P1,P2); % from projection matrices, ground truth
Fg = Fg/norm(Fg);
end