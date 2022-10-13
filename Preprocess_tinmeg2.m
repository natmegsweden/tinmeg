
%% Specify conditions, event triggers and subject list
run('cond_trig_sub_tinmeg2.m');

meg_data_path = '/archive/20061_tinnitus/MEG/';

%% Find tinmeg2-files for subject and create cell-array of file paths

% Create cell array for subjects filepaths
subpaths = cell(1);

for i = 1:height(sub_date);

% Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
fnames = find_files(subpath, {'tinmeg2', 'tsss'}, 'ds');

subpaths{i,1} = char(sub_date.ID{i}); %Include ID for tracking

    for fileindex = 1:length(fnames);
        subpaths{i,1+fileindex} = [subpath char(fnames(fileindex))]; % NB! n of files differ between rows, some subjects have empty columns
    end
clear fnames subpath i fileindex
end

%writetable(cell2table(subpaths), '../Analysis Output/included_filepaths_tinmeg2.csv') %Write log

%% create log for n of trials in raw data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_raw.csv', 'file');
    rawcondlog = readtable('../Analysis Output/n_cond_raw.csv', 'ReadVariableNames', false);
    rawcondlog = table2cell(rawcondlog);
else 
    %NB 22 columns hardcoded
    rawcondlog = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

%% Preprocess - loop for subject list in Full analysis

% specify i for now
fs_ds = 200;

%Check if ID is in TRIAL-log and determine row in log to write to
if find(strcmp(['ID' char(subpaths(i,1))], rawcondlog)) > 0;
    rawlogheight = find(strcmp(['ID' char(subpaths(i,1))], rawcondlog));
else
    %Find height of trial-log and +1 for new row
    rawlogheight = size(rawcondlog);
    rawlogheight = rawlogheight(1) + 1;

    %Write ID to new row, column 1
    rawcondlog{rawlogheight,1} = ['ID' char(subpaths(i,1))];
end

outdir = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/']



%Check if subject dir exist, create/define
if ~exist(outdir, 'file');
mkdir(outdir);
end

%For each condition
for ii = 1:length(conditions)

%Create filename and check if raw file for condition exist
fname = [conditions{ii} '_raw' '.mat'];

    if exist([outdir fname], 'file')
    warning(['Output' fname ' exist for subject ' sub_date.ID{i}])
    continue
    end

% Define trials
cfg = [];

sub_data = subpaths(i,2:max(find(~cellfun(@isempty,subpaths(i,:)))));

% NB! cellfun for cfg.dataset defines 2:highest populated column
cfg.dataset             = sub_data;
cfg.trialdef.prestim    = 0.500;        % seconds before trigger
cfg.trialdef.poststim   = 0.500;        % seconds after trigger
cfg.trialdef.eventvalue = eval(['cond.' char(conditions(ii)) 'trig']); % :/
cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';

% On pre/post-stim, consider ISI240 as longest trial, i.e trial structure of
% baseline: 200 ms
% gap duration: 50 ms
% ISI: 240 ms, i.e, min of: 490 ms prestim

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

pre_tinmeg2 = ft_preprocessing(cfg);

save([outdir fname], 'pre_tinmeg2', '-v7.3')

%write n of trials log        
%How many stimtypes for cond and what trigger values
nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
trigs = eval(['cond.' char(conditions(ii)) 'trig']);

    %omg this gets ugly..
    %Write stim to rawcondlog cell-array
    for iii = 1:length(eval(['cond.' char(conditions(ii)) 'trig']))
    if strcmp(conditions(ii), 'PO60')
    rawcondlog(rawlogheight,iii+1) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'PO70')
    rawcondlog(rawlogheight,iii+7) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP60')
    rawcondlog(rawlogheight,iii+12) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP70')
    rawcondlog(rawlogheight,iii+16) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GO')
    rawcondlog(rawlogheight,iii+20) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    end
    end

%Write log to csv
writetable(cell2table(rawcondlog), '../Analysis Output/n_cond_raw_tinmeg2.csv', 'WriteVariableNames', false) %Write log

%downsample and save
cfg = [];
cfg.resamplefs = fs_ds;

pre_tinmeg2_ds = ft_resampledata(cfg, pre_tinmeg2);

%NB new fname
fname = [char(conditions(ii)) '_ds' '.mat'];

save([outdir fname], 'pre_tinmeg2_ds')

%clear temp variables
clear pre_tinmeg2 pre_tinmeg2_ds

%end for conditions
end
