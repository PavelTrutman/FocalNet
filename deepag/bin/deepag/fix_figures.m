% fixes jpg figures into your current matlab version
% usefull when you only have 2012 matlab on cluster
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

stdir=cd;
cd bougnoux;
figs=dir('*.fig');
for i=1:size(figs,1)
    fig=openfig(figs(i).name,'invisible');
    saveas(fig,[figs(i).name(1:end-3) 'jpg']);
end
cd sdtir;