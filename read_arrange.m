%% Read and arrange fif-files for TinMEG1
%
%
%

%% Read subjects and dates

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

n_subjects = height(sub_date);
disp(['Number of subjects in table is ' num2str(n_subjects)])

%% Find relevant files for subject and create cell-array of file paths

% Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{1}) '/' char(sub_date.date{1}) '/'];
fnames = find_files(subpath, {'tinmeg1', 'tsss'});

%is exc_str for find-files working?
%order matters?

% Create cell array for subjects files
fpaths = cell(1, length(fnames));

for fileindex = 1:length(fnames);
    
    fpaths{1,fileindex} = [subpath char(fnames(fileindex))];
    
end

% Implement write to log file here i.e. n of files for subject

%% Specify event triggers

%Trigger value (STI101) at pulse onset
events60PO = [40968 36872 34824 33800 33288 33032];
events70PO = [36876 34828 33804 33292 33036];
events60GP = [49800 49736 49704 49688];
events70GP = [49804 49740 49708 49692];

%Trigger value (STI101) at gap onset
eventsGO   = [16386 16390];

% Pulse trigger:
% B_C60_P70: 40968
% B_C60_P75: 36872
% B_C60_P80: 34824
% B_C60_P85: 33800
% B_C60_P90: 33288
% B_C60_P95: 33032

% B_C70_P75: 36876
% B_C70_P80: 34828
% B_C70_P85: 33804
% B_C70_P90: 33292
% B_C70_P95: 33036

%ISI pulse trigger (gaptrig = pulsetrig - 6)
% A_C60_i0: 49800
% A_C60_i60: 49736
% A_C60_i120: 49704
% A_C60_i240: 49688

% A_70_i0: 49804
% A_70_i60: 49740
% A_70_i120: 49708
% A_70_i240: 49692

%GO gap-trig:
% A_GO_60: 16386
% A_GO_70: 16390

%% Define trials and preprocess

% Define trials
cfg = [];
cfg.dataset             = fpaths;      % Cell array of subjects files
cfg.trialdef.prestim    = 0.35;        % seconds before trigger
cfg.trialdef.poststim   = 0.30;        % seconds after trigger
cfg.trialdef.eventvalue = events60PO;
cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';

cfg = ft_definetrial(cfg);

%Remove trials from cfg.trl that have negative sample index for trial start
cfg.trl = cfg.trl(cfg.trl(:,1) >= 0,:);
  
%Remove trials from cfg.trl that have higher sample index than exist in file
cfg.trl = cfg.trl(cfg.trl(:,2) < max([cfg.event.sample]),:);

% preprocessing
cfg.demean     = 'no';
cfg.lpfilter   = 'no';
cfg.hpfilter   = 'no';
cfg.dftfilter  = 'no';
cfg.channel    = {'MEG', 'ECG', 'EOG'};

epochs = ft_preprocessing(cfg);

%expect 50
for i = 1:length(events)
    disp(events(i))
    tottrig(i) = sum(epochs.trialinfo == events(i))
end
