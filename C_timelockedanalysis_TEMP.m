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
    % I added this little snippet to control input/output file management
    infile = find_files(datapath, {sub_date.ID{i} conditions(ii), 'clean'})
    % Alternativly /find_files is meant to find raw files among many
    % possible. At later stages in the processing, the filenames hsould be
    % fixed)
    sub_dir = fullfile(datapath, sub_date.ID{i})
    infile = fullfile(sub_dir, [conditions(ii), 'clean']) % fix the strings in the [] to match actial filenames
        
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
        outfname = fullfile(sub_dir, [conditions(ii), '-epo.mat']) % fix the strings in the [] to match actial filenames
        
        trigger = trigs(stim_index);

        cfg = [];
        cfg.covariance              = 'yes';
        cfg.covariancewindow        = 'prestim';
        cfg.keeptrials              = 'no'; %if yes, no avg in output variable "timelockeds"?
        cfg.preproc.demean          = 'yes';
        cfg.preproc.baselinewindow  = [-0.200 -0.100]; %check this (!! this is incompatible with the trialdefinition from line 47)
        cfg.preproc.lpfilter        = 'yes';
        cfg.preproc.lpfreq          = 70;
        cfg.trials = temptrials.trialinfo == trigger;

        epo = ft_timelockanalysis(cfg, temptrials);
         
        % save and arrange - since this will be big struct you might
        % consider saving individually and then arrange at a later stage
        % before grand average and/or statistics.
        save(outfname, epo)    % Outname defined above
        
        epochs.(conditions{ii}){i, stim_index} = epo;

    end
    
    %end for condition
    end
%end for subject
end

clear('trigs', 'trigger', 'i', 'ii', 'nstim', 'stim_index', 'fname', 'cleaned4mat', 'temptrials');

%Large file size, save with caution
%save('../mat_data/epochs.mat', 'epochs', '-v7.3')

% Alternative: rather than saving all in one struct, consider making one
% struct for each condition. Then you only need to load the condition that
% you use for a given operation. Tis is example code of one way to do it,
% you will have to tweak it to make it work for your data.

allTaskA = []  % cells instead of struct
allTaskB = []
allTaskC = []
for ii = 1:20
    load(fullfile(pathtodata, 'dataAepo.mat')
    allTaskA{ii} = epo
    clear epo % just in case the object has the same name
    
    load(fullfile(pathtodata, 'dataBepo.mat')
    allTaskB{ii} = epo
    
    load(fullfile(pathtodata, 'dataBepo.mat')
    allTaskC{ii} = epo
    
    % and so on...
end

save(fullfile(datapath, 'allTaskA.mat'), allTaskA)
save(fullfile(datapath, 'allTaskB.mat'), allTaskB)
save(fullfile(datapath, 'allTaskC.mat'), allTaskC)

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
