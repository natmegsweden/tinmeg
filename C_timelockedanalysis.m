%%

% To do:
% Check if exist and skip
% Plot

%% Find and load variables and files

load('../mat_data/conditions.mat');
load('../mat_data/cond.mat');

%Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

disp(['Number of subjects in table is ' num2str(height(sub_date))])

%Only one condition here so far
datapath = '../mat_data/';

%Create structure to write ouput

epochs = struct;

%% Create timelockeds

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

%find ID:
%find(contains(epochs.ID, 'texthere'))
