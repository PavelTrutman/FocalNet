% A script for generating scatter plot series of Bougnoux formula results
% 
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

% 

plotting=true;
calculating=true;
deepagpaths;

if calculating
    % set parameters
    file='../../data/paris/correspondences_F_synth_1K.mat';
    %file='../../data/paris/correspondences_F.mat';
    corr=[7:17 18:3:29 30:5:49 50:10:100];  % numbers of correspondences - x axis
    corr=[8:12];
    %corr=[8];
    pop_size=100; % number of points that will be on scatter plot
    n=size(corr,2); % x axis length
    stat_data=NaN(4,n,4);
        
    % size(stat_data):=(methods, n, parameters), where:
    % methods are as in the legend below
    % parameters are [absolute_error a_std relative_error r_std]
    for i=1:n
        Fparam.corr=corr(i);
        if corr(i)>7
            % beware, data are magically reshaped here
            stat_data(1:2,i,:)=permute(bougnoux_scatter(file,corr(i),pop_size,'Free'),[1 3 2]);
        end
        stat_data(3:4,i,:)=permute(bougnoux_scatter(file,corr(i),pop_size,'|F|=0'),[1 3 2]);
    end
end

if plotting
    titles(1)={'mean error of estimating (f1,f2)'};
    legends(1,:)={'only real f', 'all f','singularized F, only real f', 'singularized F, all f'};
    ylabels(1)={'error, px'};
    names(1)={'bougnoux/error'};
    
    titles(2)={'std of estimating (f1,f2)'};
    legends(2,:)=legends(1,:);
    ylabels(2)={'std, px'};
    names(2)={'bougnoux/std'};
    
    titles(3)={'mean error of estimating f1/f2'};
    legends(3,:)=legends(1,:);
    ylabels(3)={'error'};
    names(3)={[names{1} '_proportion']};
    
    titles(4)={'std of estimating f1/f2'};
    legends(4,:)=legends(1,:);
    ylabels(4)={'std'};
    names(4)={[names{2} '_proportion']};
    
    for i=1:4
        plot(corr,abs(stat_data(1,:,i)), 'bx-',corr,abs(stat_data(2,:,i)),'rx-',corr,abs(stat_data(3,:,i)), 'mx-',corr,abs(stat_data(4,:,i)),'kx-');
        legend(legends(i,:));
        title(titles(i));
        xlabel('number of correspondences');
        ylabel(ylabels(i));
        saveas(gcf,[names{i} '.fig']);
        saveas(gcf,[names{i} '.jpg']);
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperPosition', [0 0 1024 800]);
        set(gcf,'Position',[100, 100, 1024, 800]);
        hold off
    end
end


