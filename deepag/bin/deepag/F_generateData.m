% Script for generating synthetic Fundamental matrix feature vectors
%
% saves in corresponding files simulation of F and u. Files then could be
% used in scripts bougnoux_scatter, k_nn, nn_*
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016
function F_generateData()
per_corr=1; % samples (with different noise) per correspondence set
corrnum=20; % number of correspondences in simulation
tic();
trSize = 100*1000;
valSize = 1*100;
train = 0; % whether the set should be training or validating
saveFeatures=false;
saveCorrs=true;
minlen=200; diflen=2000; % focal len
noise=1; % noise in correspondences

if train
    Size=trSize;
else
    Size=valSize;
end

coefs=zeros(9,Size*per_corr);
f=zeros(2,Size*per_corr);
norm_=zeros(2,Size*per_corr);
U=zeros(4*corrnum,Size*per_corr);
for i=1:Size
    % generate random features, reshape it (to vector) 
    [F,f1,f2,A,u]=getRandomF(minlen,diflen,noise,per_corr,corrnum);
    for j=1:size(F,1)
        index=(i-1)*per_corr+j;
        repS = adprintf({}, [num2str(index), '/', num2str(Size*per_corr)]);
        coefs(:,index)=reshape(F{j},9,1);
        norm_(:,index)=[1/A{j}{1}(1) 1/A{j}{2}(1)];
        f(:,index)=[f1; f2]./norm_(:,index);
        rmprintf(repS);
        U(:,index) = [reshape(u{1}, [], 1); reshape(u{2}, [], 1)];        
    end
end

% saving
if saveFeatures
    if train
        tr_coefs=coefs;
        tr_norm=norm_;
        tr_f=f;
        save('../../data/paris/features_F_synth.mat', 'tr_coefs', 'tr_norm', 'tr_f', '-v7.3');
        clear tr_f tr_coefs tr_norm;
    else
        val_coefs=coefs;
        val_norm=norm_;
        val_f=f;
        save('../../data/paris/features_F_nsynthrep_sample_1K.mat', 'val_coefs', 'val_norm', 'val_f', '-v7.3');
        clear val_f val_coefs val_norm;
    end
end
if saveCorrs
    size_=size(f,2);
    size_3=floor(size_/3);
    ord=[0,size_3,size_3*2,size_];
    corr_tr=form_corrstructure(f,U,norm_,ord(1:2));
    corr_val=form_corrstructure(f,U,norm_,ord(2:3));
    corr_tst=form_corrstructure(f,U,norm_,ord(3:4));
    save('../../data/paris/correspondences_F_synth_1K.mat', 'corr_tr', 'corr_val', 'corr_tst', '-v7.3');
    clear corr_tr corr_val corr_tst;
end

toc();
end

function str=form_corrstructure(f,u,norm_,part)
% form structure in the convential format of stored data

str.f=f(:,(part(1)+1):part(2));
str.u=u(:,(part(1)+1):part(2));
str.norm=norm_(:,(part(1)+1):part(2));

end

function [F,f1,f2,A, u]=getRandomF(minlen,diflen,noise,per_corr,corrnum)
%get random cell array of matrices F 
% F - fundamental matrix
% f1,f2 - focal lengths
% A - scaling matrices
% u - correspondences, beware, has different format

Fparam.f1 = minlen+diflen*rand(1,1);
Fparam.f2 = minlen+diflen*rand(1,1);
Fparam.per_corr=per_corr;
Fparam.noise=noise;
Fparam.corr=corrnum;

[F,A,u]=F_simulate(Fparam,'Free');
f1=Fparam.f1;
f2=Fparam.f2;
end
