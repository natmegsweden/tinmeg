
%Load conditions for subjects, average trials via ft_timelockanalysis

sensspace_avg = struct;

for i = 1:4%length(sub_date.ID)
    
    subinpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    sensspace_avg.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1%:length(conditions)
    
    tempdat = load([subinpath conditions{ii} 'ica.mat']);
    tempdat = tempdat.([conditions{ii} 'ica']);
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    trig = eval(['cond.' char(conditions(ii)) 'trig']);

    
        for iii = 1:nstim

        cfg = [];
        cfg.covariance = 'yes';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no' %if yes, no avg in output variable "timelockeds"?
        cfg.preproc.demean = 'yes';
        cfg.preproc.baselinewindow = [-0.200 -0.100]; %check this
        cfg.preproc.lpfilter = 'yes';
        cfg.preproc.lpfreq = 70;
        cfg.trials = tempdat.trialinfo == trig(iii);

        sensspace_avg.(conditions{ii}){i, iii} = ft_timelockanalysis(cfg, tempdat);
        
        %For trigs
        end
    
    %For conditions
    end

%For subjects
end



grandavg = struct();

for i = 1:2%numel(sub_date.ID)
   
    grandavg(i).PO60 = sensspace_avg.PO60{i,5}.avg;
    
end

mean_sub = sensspace_avg.PO60{1,1};
mean_sub.avg = mean(cat(3, grandavg.PO60), 3);

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.linecolor = 'brgkym'
ft_multiplotER(cfg, sensspace_avg.PO60{1:2,5}, mean_sub);


