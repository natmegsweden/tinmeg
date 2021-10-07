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
% if find(strcmp(['ID' char(subpaths(i,1))], artcondlog)) > 0;
%     artlogheight = find(strcmp(['ID' char(subpaths(i,1))], artcondlog));
% else
%     %Find height of trial-log and +1 for new row
%     artlogheight = size(artcondlog);
%     artlogheight = artlogheight(1) + 1;
% 
%     %Write ID to new row, column 1
%     artcondlog{artlogheight,1} = ['ID' char(subpaths(i,1))];
% end

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
cfg.trialdef.poststim   = 0.320;        % seconds after trigger
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

%offset PO-trials with 50ms to compensate for ERROR in presentation-script that had trigger wait for gap duration (50ms) even when no gap
for r = 1:length(cfg.trl)
    if (ismember(cfg.trl(r,4), [cond.PO70trig cond.PO60trig]))
        cfg.trl(r,1) = cfg.trl(r,1) - 250;
        cfg.trl(r,2) = cfg.trl(r,2) - 250;
        cfg.trl(r,3) = cfg.trl(r,3) - 250;
    end
end

clear r

% preprocessing
cfg.demean     = 'no';
cfg.lpfilter   = 'no';
cfg.hpfilter   = 'no';
cfg.dftfilter  = 'no';
cfg.channel    = {'MEG', 'ECG', 'EOG'};

res4mat = ft_preprocessing(cfg);

%adjust time-variable to match offset for PO-trials
if conditions{ii} == 'PO60' | conditions{ii} == 'PO70'    
    cfg = [];
    cfg.offset = 250;
    res4mat = ft_redefinetrial(cfg, res4mat);   
end

save([outdir fname], 'res4mat', '-v7.3')

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

%NB new fname
fname = [char(conditions(ii)) '_ds' '.mat'];

save([outdir fname], 'res4mat_ds')

%clear temp variables
clear res4mat res4mat_ds

%end for conditions
end
