
%Load conditions for subjects, average trials via ft_timelockanalysis
sensspace_avg = struct;
all_avg = struct();
gravg = struct();
for i = 1:4%length(sub_date.ID)
    
    subinpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    sensspace_avg.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    tempdat = load([subinpath conditions{ii} 'ica.mat']);
    tempdat = tempdat.([conditions{ii} 'ica']);
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    trig = eval(['cond.' char(conditions(ii)) 'trig']);

    
        for iii = 1:nstim

        cfg = [];
        cfg.covariance = 'yes';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"?
        cfg.preproc.demean = 'yes';
        cfg.preproc.baselinewindow = [-0.200 -0.100]; %check this
        cfg.preproc.lpfilter = 'yes';
        cfg.preproc.lpfreq = 70;
        cfg.trials = tempdat.trialinfo == trig(iii);

        sensspace_avg.(conditions{ii}){i, iii} = ft_timelockanalysis(cfg, tempdat);
        
        all_avg.(conditions{ii}){i, iii} = sensspace_avg.(conditions{ii}){i, iii}.avg;
        
        gravg.(conditions{ii}){iii} = mean(cat(3, all_avg.(conditions{ii}){:, iii}), 3);
        
        %For trigs
        end
    
    %For conditions
    end

%For subjects
end

save('../mat_data/sensspace/sensspace_avg.mat', 'sensspace_avg', '-v7.3');
save('../mat_data/sensspace/all_avg.mat', 'all_avg', '-v7.3');
save('../mat_data/sensspace/grand_avg.mat', 'gravg', '-v7.3');

%% MultiplotER
mean_sub = sensspace_avg.PO60{1,1};
mean_sub.avg = gravg.PO60{5};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.linecolor = [repmat(.5, 4, 3); 1 0 0]; %change according to n of subjects
ft_multiplotER(cfg, sensspace_avg.PO60{1:4,5}, mean_sub);

%% TopoplotER
mean_sub = sensspace_avg.PO60{1,1};

    cfg = [];
    cfg.parameter = 'avg';
    cfg.layout = 'neuromag306mag.lay';
    cfg.colorbar = 'no';
    cfg.comment = '';

    figure;
    
for i = 1:numel(cond.PO60label);
    
    mean_sub.avg = gravg.PO60{i};

    subplot(2,3,i);
    ft_topoplotER(cfg, mean_sub);
    %title(cond.PO60label{i}, 'Interpreter', 'none');
    
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

end
