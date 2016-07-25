% deepag.m - Deep Learning 4 Algebraic Geometry
% T. Pajdla, pajdla@gmail.cz
% 1 Jan 2016 - 31 Dec 2016

%close all
drawnow
%% Base path
pw = pwd; cd(fullfile('..','..')); pc.basePath = pwd; cd(pw); clear('pw');
%%
pc.dataPath = fullfile(pc.basePath,'data');
run(fullfile(pc.dataPath,'deepagparini.m')); % default parameters
% data specific paramaters
fn = {fullfile(pc.dataPath,ps.Data,'deepagpar.m')};
for i=1:numel(fn)
    if exist(fn{i},'file')
        fid=fopen(fn{i},'r');
        while true
            s=fgets(fid);
            if ~ischar(s),
                break;
            else
                eval(s);
            end
        end
        fclose(fid);
    end
end
%
pc.clrs = ['k','b','g','r','m','c']; % colors
pc.lsts = {'.-','s-','.:','.-.','-','--'}; % line styles
%
clear fid fn i s pw 

