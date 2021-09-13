%%

% To do:
% Check if exist and skip
% Plot

%% Create timelockeds for EOG

%implement: if exist load

for i = 1:length(sub_date.ID)

%for each condition
for ii = 1:length(conditions)

%NB reading in EOG data from file BEFORE ft_rejectvisual
fname = [conditions{ii} '_ds' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

%loads as 'res4mat_ds'
load(fpath);

%Different TOI for PO and GP compensating 50ms (gap length) for PO-trials
if ismember(conditions{ii}, {'PO60', 'PO70'})
    toilow = -0.150
    toihigh = 0.250
elseif ismember(conditions{ii}, {'GP60', 'GP70', 'GO'})
    toilow = -0.100
    toihigh = 0.300
end

%Define toi for trial in file
cfg = [];
cfg.toilim = [toilow toihigh];

temptrials = ft_redefinetrial(cfg, res4mat_ds);

%Find triggers and stims in condition
nstim = length(eval(['cond.' conditions{ii} 'trig']));
trigs = eval(['cond.' conditions{ii} 'trig']);

%Extract EOG-timelockeds for all stim types
    for stim_index = 1:nstim

    trigger = trigs(stim_index);
    label = eval(['cond.' conditions{ii} 'label']);

    cfg = [];

    cfg.channel = {'EOG001', 'EOG002'};

    cfg.covariance = 'yes';
    cfg.covariancewindow = 'prestim';
    cfg.keeptrials = 'no' %if yes, no avg in output variable "timelockeds"?
    cfg.preproc.demean = 'yes';
    cfg.preproc.baselinewindow = [toilow-0.100 toilow]; %check this
    cfg.preproc.lpfilter = 'yes';
    cfg.preproc.lpfreq = 70;
    cfg.trials = temptrials.trialinfo == trigger;
    
    eog_timelockeds = ft_timelockanalysis(cfg, temptrials);
    
    destdirectory = ['../mat_data/timelockeds/' 'ID' sub_date.ID{i} '/'];
    
    if ~exist(destdirectory, 'file');
    mkdir(destdirectory);
    
    %??
    save([destdirectory label{stim_index} '_eog' '.mat'], 'eog_timelockeds');
    
    continue
    end
    
    save([destdirectory label{stim_index} '_eog' '.mat'], 'eog_timelockeds');
    
    %end all stimtypes
    end

%end all conditions
end

%end for subject in sub_date
end


%% EOG, Keep trials = yes, (only EOG002)

for i = 1:length(sub_date.ID)

%for each condition
for ii = 1:length(conditions)

%NB reading in EOG data from file BEFORE ft_rejectvisual
fname = [conditions{ii} '_ds' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

%loads as 'res4mat_ds'
load(fpath);

%Different TOI for PO and GP compensating 50ms (gap length) for PO-trials
if ismember(conditions{ii}, {'PO60', 'PO70'})
    toilow = -0.150
    toihigh = 0.250
elseif ismember(conditions{ii}, {'GP60', 'GP70', 'GO'})
    toilow = -0.100
    toihigh = 0.300
end

%Define toi for trial in file
cfg = [];
cfg.toilim = [toilow toihigh];

temptrials = ft_redefinetrial(cfg, res4mat_ds);

%Find triggers and stims in condition
nstim = length(eval(['cond.' conditions{ii} 'trig']));
trigs = eval(['cond.' conditions{ii} 'trig']);

%Extract EOG-timelockeds for all stim types
    for stim_index = 1:nstim

    trigger = trigs(stim_index);
    label = eval(['cond.' conditions{ii} 'label']);

    cfg = [];

    cfg.channel = {'EOG002'};

    cfg.covariance = 'yes';
    cfg.covariancewindow = 'prestim';
    cfg.keeptrials = 'yes' %if yes, no avg in output variable "timelockeds"?
    cfg.preproc.demean = 'yes';
    cfg.preproc.baselinewindow = [toilow-0.100 toilow]; %check this
    cfg.preproc.lpfilter = 'yes';
    cfg.preproc.lpfreq = 70;
    cfg.trials = temptrials.trialinfo == trigger;
    
    eog_timelockeds_all = ft_timelockanalysis(cfg, temptrials);
    
    destdirectory = ['../mat_data/timelockeds/' 'ID' sub_date.ID{i} '/'];
    
    if ~exist(destdirectory, 'file');
    mkdir(destdirectory);
    
    %??
    save([destdirectory label{stim_index} '_eog_all' '.mat'], 'eog_timelockeds_all');
    
    continue
    end
    
    save([destdirectory label{stim_index} '_eog_all' '.mat'], 'eog_timelockeds_all');
    
    %end all stimtypes
    end

%end all conditions
end

%end for subject in sub_date
end

%%

figure
hold on
plot(epochs_eog.PO60{1,1}.avg(2,:))
plot(epochs_eog.PO60{1,2}.avg(2,:))
plot(epochs_eog.PO60{1,3}.avg(2,:))
plot(epochs_eog.PO60{1,4}.avg(2,:))
plot(epochs_eog.PO60{1,5}.avg(2,:))
plot(epochs_eog.PO60{1,6}.avg(2,:))
hold off

