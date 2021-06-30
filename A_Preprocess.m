%% Read and arrange fif-files for TinMEG1
%
% To do or consider:
% Append filenames to included_filepaths.csv - currently overwrites (include date?)
% Modify cfg.trialdef.pre/poststim for PO/GP trials

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

%Check if ID is in ARTIFACT-log and determine row in log to write to
if find(strcmp(['ID' char(subpaths(i,1))], artcondlog)) > 0;
    artlogheight = find(strcmp(['ID' char(subpaths(i,1))], artcondlog));
else
    %Find height of trial-log and +1 for new row
    artlogheight = size(artcondlog);
    artlogheight = artlogheight(1) + 1;

    %Write ID to new row, column 1
    artcondlog{artlogheight,1} = ['ID' char(subpaths(i,1))];
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

sub_data = subpaths(i,2:max(find(~cellfun(@isempty,subpaths(i,:)))));

% NB! cellfun for cfg.dataset defines 2:highest populated column
cfg.dataset             = sub_data;
cfg.trialdef.prestim    = 0.35;        % seconds before trigger
cfg.trialdef.poststim   = 0.30;        % seconds after trigger
cfg.trialdef.eventvalue = eval(['cond.' char(conditions(ii)) 'trig']); % :/
cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';

cfg = ft_definetrial(cfg);

%Remove trials from cfg.trl that have negative sample index for trial start
cfg.trl = cfg.trl(cfg.trl(:,1) >= 0,:);

%Remove trials from cfg.trl that have higher sample index than exist in file
cfg.trl = cfg.trl(cfg.trl(:,2) < max([cfg.event.sample]),:);



%Removed artifact script here

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
writetable(cell2table(rawcondlog), '../Analysis Output/n_cond_raw.csv', 'WriteVariableNames', false) %Write log

%downsample and save
cfg = [];
cfg.resamplefs = fs_ds;

res4mat_ds = ft_resampledata(cfg, res4mat);

fname = ['ID' char(subpaths(i,1)) '_' char(conditions(ii)) '_ds' '.mat'];
fpath = ['../mat_data/' fname];

save(fpath, 'res4mat_ds')

%clear temp variables
clear res4mat res4mat_ds
end
