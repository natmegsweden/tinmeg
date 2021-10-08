%% Create timelockeds for EOG

for i = 1%:length(sub_date.ID)

destdirectory = ['../mat_data/timelockeds/' 'ID' sub_date.ID{i} '/EOG/'];

%Check if output folder exist and skip subject or create folder
if exist(destdirectory, 'file');
   warning(['EOG output folder already exist for subject ' sub_date.ID{i} ' - skipping..']);
    continue
elseif ~exist(destdirectory, 'file');
    mkdir(destdirectory);
end
    
    %for each condition
    for ii = 1:length(conditions)

    %NB reading in EOG data from file BEFORE ft_rejectvisual
    fname_in = [conditions{ii} '_ds' '.mat'];
    fpath_in = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname_in];

    %loads as 'res4mat_ds'
    load(fpath_in);

    %Find triggers and stims in condition
    nstim = length(eval(['cond.' conditions{ii} 'trig']));
    trigs = eval(['cond.' conditions{ii} 'trig']);

    %Extract EOG-timelockeds for all stim types
        for stim_index = 1:nstim

        trigger = trigs(stim_index);
        label = eval(['cond.' conditions{ii} 'label']);

        fname_out = [label{stim_index} '_eog' '.mat'];

        cfg = [];

        cfg.channel = {'EOG002'};

        cfg.covariance = 'yes';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"?
        cfg.preproc.demean = 'yes';

        if ismember(conditions{ii}, {'GO', 'PO60', 'PO70'});
        %Baseline window: 150ms before pulse onset in PO trials (gap onset for GO)
        cfg.preproc.baselinewindow = [-0.150 0];
        elseif ismember(conditions{ii}, {'GP60', 'GP70'});
            %Baseline window variable timepoint: 150ms before gap onset in GP trials
            if stim_index == 1 %ISI 0
            cfg.preproc.baselinewindow = [-0.200 -0.050];
            elseif stim_index == 2 %ISI 60
                cfg.preproc.baselinewindow = [-0.260 -0.110];
            elseif stim_index == 3 %ISI 120
                cfg.preproc.baselinewindow = [-0.320 -0.170];
            elseif stim_index == 4 %ISI 240
                cfg.preproc.baselinewindow = [-0.440 -0.290];
            end
        end

        cfg.preproc.lpfilter = 'yes';
        cfg.preproc.lpfreq = 70;
        cfg.trials = res4mat_ds.trialinfo == trigger;

        eog_timelockeds = ft_timelockanalysis(cfg, res4mat_ds);

        save([destdirectory fname_out], 'eog_timelockeds');
        
        clear eog_timelockeds
        
        %run again with cfg.keeptrials = yes
        cfg.keeptrials = 'yes';
        eog_timelockeds_all = ft_timelockanalysis(cfg, res4mat_ds);
        
        fname_out = [label{stim_index} '_eog_all' '.mat'];
        
        save([destdirectory fname_out], 'eog_timelockeds_all');
        
        clear eog_timelockeds_all
        
        %end all stimtypes
        end

    %end all conditions
    end

%end for subject in sub_date
end
