%% Read and arrange fif-files for TinMEG1
%
% To do or consider:
% Append filenames to included_filepaths.csv - currently overwrites (include date?)
% Modify cfg.trialdef.pre/poststim for PO/GP trials

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

trl = cfg.trl;

%Artifact rejection!
cfg = [];
cfg.trl = trl;
cfg.datafile = subpaths(i,2:max(find(~cellfun(@isempty,subpaths(i,:)))));
cfg.headerfile = subpaths(i,2:max(find(~cellfun(@isempty,subpaths(i,:)))));
cfg.continuous = 'yes';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = 'MEG';
cfg.artfctdef.zvalue.cutoff = 20;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.cumulative = 'yes';
cfg.artfctdef.zvalue.medianfilter = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff = 'yes';

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_jump] = ft_artifact_zvalue(cfg);


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
cfg.resamplefs = fs_ds;

res4mat_ds = ft_resampledata(cfg, res4mat);

fname = ['ID' char(subpaths(i,1)) '_' char(conditions(ii)) '_ds' '.mat'];
fpath = ['../mat_data/' fname];

save(fpath, 'res4mat_ds')

%clear temp variables
clear res4mat res4mat_ds
end
