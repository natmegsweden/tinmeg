%% Read and arrange fif-files for TinMEG1
%
% To do or consider:
% Append filenames to included_filepaths.csv - currently overwrites (include date?)
% Modify cfg.trialdef.pre/poststim for PO/GP trials
%
%

%% Read subjects and dates

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

disp(['Number of subjects in table is ' num2str(height(sub_date))])

%% Find relevant files for subject and create cell-array of file paths

% Create cell array for subjects filepaths
subpaths = cell(1);

for i = 1:height(sub_date);

% Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
fnames = find_files(subpath, {'tinmeg1', 'tsss'}, 'ds');

subpaths{i,1} = char(sub_date.ID{i}); %Include ID for tracking

    for fileindex = 1:length(fnames);
        subpaths{i,1+fileindex} = [subpath char(fnames(fileindex))]; % NB! n of files differ between rows, some subjects have empty columns
    end
clear fnames subpath i fileindex
end

writetable(cell2table(subpaths), '../Analysis Output/included_filepaths.csv') %Write log

%% Specify conditions and event triggers

%all conditions
conditions = ({'PO60', 'PO70', 'GP60', 'GP70', 'GO'});
%Save structure of conditions
save('../mat_data/conditions.mat', 'conditions');

%Structure for triggers and labels
cond = struct();

%trigger at pulse onset
cond.PO60trig   = [40968 36872 34824 33800 33288 33032];
cond.PO60label  = ({'PO60_70', 'PO60_75', 'PO60_80', 'PO60_85', 'PO60_90', 'PO60_95'});
cond.PO70trig   = [36876 34828 33804 33292 33036];
cond.PO70label  = ({'PO70_75', 'PO70_80', 'PO70_85', 'PO70_90', 'PO70_95'});
cond.GP60trig   = [49800 49736 49704 49688];
cond.GP60label  = ({'GP60_i0', 'GP60_i60', 'GP60_i120', 'GP60_i240'});
cond.GP70trig   = [49804 49740 49708 49692];
cond.GP70label  = ({'GP70_i0', 'GP70_i60', 'GP70_i120', 'GP70_i240'});

%trigger at gap onset (gap only)
cond.GOtrig     = [16386 16390];
cond.GOlabel    = ({'GO_60', 'GO_70'});

save('../mat_data/cond.mat', 'cond');

%ISI pulse trigger (gaptrig = pulsetrig - 6)
% A_C60_i0: 49800
% A_C60_i60: 49736
% A_C60_i120: 49704
% A_C60_i240: 49688

% A_70_i0: 49804
% A_70_i60: 49740
% A_70_i120: 49708
% A_70_i240: 49692

%% Define trials and preprocess conditions

%create log for n of trials in raw data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_raw.csv', 'file');
    rawcondlog = readtable('../Analysis Output/n_cond_raw.csv', 'ReadVariableNames', false);
    rawcondlog = table2cell(rawcondlog);
else 
    %NB 22 columns hardcoded
    rawcondlogl = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

%For each subject in sub_date table
for i = 1:height(sub_date)
    
    %Check if ID is in trial-log and determine row in log to write to
    if find(strcmp(['ID' char(subpaths(i,1))], rawcondlog)) > 0;
        logheight = find(strcmp(['ID' char(subpaths(i,1))], rawcondlog));
    else
        %Find height of trial-log and +1 for new row
        logheight = size(rawcondlog);
        logheight = logheight(1) + 1;
  
        %Write ID to new row, column 1
        rawcondlog{logheight,1} = ['ID' char(subpaths(i,1))];
    end
    
%For each event in condevents
for ii = 1:length(conditions)
    
    %Create filename and check if raw file for condition exist
    fname = ['ID' char(subpaths(i,1)) '_' char(conditions(ii)) '_raw' '.mat'];
    fpath = ['../mat_data/' fname];
    
    if exist(fpath, 'file')
    warning(['Output' fname ' exist for subject ' char(subpaths(i,1))])
    continue
    end
    
    % Define trials
    cfg = [];

    % NB! cellfun for cfg.dataset defines 2:highest populated column
    cfg.dataset             = subpaths(i,2:max(find(~cellfun(@isempty,subpaths(i,:)))));
    cfg.trialdef.prestim    = 0.35;        % seconds before trigger
    cfg.trialdef.poststim   = 0.30;        % seconds after trigger
    cfg.trialdef.eventvalue = eval(['cond.' char(conditions(ii)) 'trig']); % :/
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

    res4mat = ft_preprocessing(cfg);

    save(fpath, 'res4mat', '-v7.3')

    %write n of trials log        
    %How many stimtypes for cond and what trigger values
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    trigs = eval(['cond.' char(conditions(ii)) 'trig']);
    
    %omg this gets ugly..
    %Write stim to rawcondlog cell-array
    for iii = 1:length(eval(['cond.' char(conditions(ii)) 'trig']))
    if strcmp(conditions(ii), 'PO60')
    rawcondlog(logheight,iii+1) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'PO70')
    rawcondlog(logheight,iii+7) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP60')
    rawcondlog(logheight,iii+12) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP70')
    rawcondlog(logheight,iii+16) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GO')
    rawcondlog(logheight,iii+20) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    end
    end
    
    %Write log to csv
    writetable(cell2table(rawcondlog), '../Analysis Output/n_cond_raw.csv', 'WriteVariableNames', false) %Write log
    
    %downsample and save
    cfg = [];
    cfg.resamplefs = 200;

    res4mat_ds = ft_resampledata(cfg, res4mat);

    fname = ['ID' char(subpaths(i,1)) '_' char(conditions(ii)) '_ds' '.mat'];
    fpath = ['../mat_data/' fname];

    save(fpath, 'res4mat_ds')

    %clear temp variables
    clear res4mat res4mat_ds

    end
    
end