%%

% To do:
% Check if exist and skip
% Plot

%% Create timelockeds for EOG

%implement: if exist load

%for each condition
for ii = 1:length(conditions)

%NB reading in EOG data from file BEFORE ft_rejectvisual
fname = ['ID' char(subpaths(i,1)) '_' char(conditions(ii)) '_ds' '.mat'];
fpath = ['../mat_data/' fname];

%loads as 'res4mat_ds'
load(fpath);

%Different TOI for PO and GP
if ismember(conditions{ii}, {'PO60', 'PO70', 'GO'})
    toilow = -0.100
    toihigh = 0.300
elseif ismember(conditions{ii}, {'GP60', 'GP70'})
    toilow = -0.100
    toihigh = 0.300
end

%Define toi for trial in file
cfg = [];
cfg.toilim = [toilow toihigh];

temptrials = ft_redefinetrial(cfg, res4mat_ds);

%Find triggers and stims in condition
nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
trigs = eval(['cond.' char(conditions(ii)) 'trig']);

%Extract EOG-timelockeds for all stim types
    for stim_index = 1:nstim

    trigger = trigs(stim_index);
    label = eval(['cond.' char(conditions(ii)) 'label']);

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
    
    destdirectory = ['../mat_data/timelockeds/' 'ID' char(subpaths(i,1)) '/'];
    
    if ~exist(destdirectory, 'file');
    mkdir(destdirectory);
    continue
    end
    
    save([destdirectory char(label(stim_index)) '_eog' '.mat'], 'eog_timelockeds');
    
    %epochs_eog.(conditions{ii}){i, stim_index} = ft_timelockanalysis(cfg, temptrials);
    
    %end all stimtypes
    end

%end all conditions
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








%For subject
for i = 1:length(sub_date.ID)
    
    %write ID to epochs output
    epochs.ID{i,1} = sub_date.ID{i};
    
    %for conditions
    for ii = 1:length(conditions);
    fname = find_files(datapath, {sub_date.ID{i} conditions(ii), 'clean'});
    
    %loads variable cleaned4mat with a *poof*
    load([char(datapath) char(fname)]);
    
    %Define toi for trial in file
    cfg = [];
    cfg.toilim = [-0.100 0.300];

    temptrials = ft_redefinetrial(cfg, cleaned4mat);
    
    %Find triggers and stims in condition
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    trigs = eval(['cond.' char(conditions(ii)) 'trig']);
    
    %Run timelockeds for all stim types
    for stim_index = 1:nstim
    
    trigger = trigs(stim_index);
    
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'prestim';
    cfg.keeptrials = 'no' %if yes, no avg in output variable "timelockeds"?
    cfg.preproc.demean = 'yes';
    cfg.preproc.baselinewindow = [-0.200 -0.100]; %check this
    cfg.preproc.lpfilter = 'yes';
    cfg.preproc.lpfreq = 70;
    cfg.trials = temptrials.trialinfo == trigger;

    epochs.(conditions{ii}){i, stim_index} = ft_timelockanalysis(cfg, temptrials);
    
    end
    
    %end for condition
    end
%end for subject
end

clear('trigs', 'trigger', 'i', 'ii', 'nstim', 'stim_index', 'fname', 'cleaned4mat', 'temptrials');

%Large file size, save with caution
%save('../mat_data/epochs.mat', 'epochs', '-v7.3')

%find ID:
%find(contains(epochs.ID, 'text_here'))

%% Plot for sanity check

% Sensor selection based on pilot only

%Timelockeds for left side MEGMAG of interest
MAGoiL = {'MEG1621', 'MEG1631', 'MEG1811', 'MEG1821', 'MEG1841'};
MAGoiR = {'MEG2221', 'MEG2411', 'MEG2421', 'MEG2431', 'MEG2441'};

%Timelockeds for right side MEGGRAD of interest
GRADoiL = {'MEG0232', 'MEG0233', 'MEG0242', 'MEG0243', 'MEG1612', 'MEG1613'};
GRADoiR = {'MEG1332', 'MEG1333', 'MEG1342', 'MEG1343', 'MEG2422', 'MEG2423'};

cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.linecolor = 'brgkym'
cfg.channel = MAGoiL
ft_singleplotER(cfg, epochs.PO60{1,1:6});
