% Script for generating synthetic Fundamental matrix feature vectors
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016
function F_features()
tic();
trSize = 100*1000;
valSize = 1*1;
train = 0; % whether the set should be training or validating
minlen=400; diflen=1000; noise=0;

if train
    tr_coefs=zeros(9,trSize);
    tr_f=zeros(2,trSize);
    tr_norm=zeros(2,trSize);
    for i=1:trSize
        [Fg,f1,f2]=getRandomF(minlen,diflen,noise);
        tr_coefs(:,i)=reshape(Fg,9,1);
        tr_f(:,i)=[f1 f2];
        tr_norm(:,i)=[1 1];
    end
    save('../../data/paris/features_F_synth.mat', 'tr_coefs', 'tr_norm', 'tr_f', '-v7.3');
    clear tr_f tr_coefs tr_norm;
else
    val_coefs=zeros(9,valSize);
    val_f=zeros(2,valSize);
    val_norm=zeros(2,valSize);
    for i=1:valSize
        [Fg,f1,f2]=getRandomF(minlen,diflen,noise);
        val_coefs(:,i)=reshape(Fg,9,1);
        f1
        f2
        F2f1f2(reshape(Fg,3,3))
        val_f(:,i)=[f1 f2];
        val_norm(:,i)=[1 1];
    end
    save('../../data/paris/features_F_nsynth_sample.mat', 'val_coefs', 'val_norm', 'val_f', '-v7.3');
    clear val_f val_coefs val_norm;
end
toc();
end

function f=getRandomF_n(P1,P2,sc,noise)
% 3D points
X = 2*sc*2*rand(3,8)+repmat([-sc;-sc;2*sc],1,8);
% Image points
u1 = X2u(X,P1);
u2 = X2u(X,P2);
% add image noise
u1(1:2,:)=u1(1:2,:)+noise*randn(size(u1(1:2,:)));
u2(1:2,:)=u2(1:2,:)+noise*randn(size(u2(1:2,:)));
% Fundamental matrix
[F,A]  = uu2F({u1,u2},{'[-1,1]','Free'}); % from image measurements
F = F/norm(F);
% Compute the difference in normalized image coordinates
Fn = inv(A{2})'*F*inv(A{1});
Fn = Fn/norm(Fn);
% Form the feature vector and the ground truth feature vector
f = Fn(:)/norm(Fn,'fro'); [~,mi] = max(abs(f)); f = f*sign(f(mi));
% Bougnoux formula aplied on normalized data and recomputed to original data
fn = F2f1f2(Fn);
fnb = fn*diag([1/A{1}(1) 1/A{2}(1)]);
end


function [F,f1,f2]=getRandomF(minlen,diflen,noise)
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
if true
    F=getRandomF_n(P1,P2,sc,noise);
    return;
else
    % Fundamental matrix
    F = PP2F(P1,P2); % from projection matrices, ground truth
    F = F/norm(F);    
end
end