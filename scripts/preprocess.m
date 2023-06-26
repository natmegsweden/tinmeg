%% Read and arrange fif-files for TinMEG1

outdir = ['../processed_data/preprocessed/' 'ID' sub_date.ID{i} '/'];

%Check if subject dir exist, create/define
if ~exist(outdir, 'file');
mkdir(outdir);
end

%For each condition
for ii = 1:length(temp_cond)

%Create filename and check if raw file for condition exist
fname = [temp_cond{ii} '_raw' '.mat'];

    if exist([outdir fname], 'file')
    warning(['Output' fname ' exist for subject ' sub_date.ID{i}])
    continue
    end

%Find what row the subjects paths are in
sub_path_row = find(strcmp(subpaths,sub_date.ID{i}));

% Define trials
cfg = [];

sub_data = subpaths(sub_path_row,2:max(find(~cellfun(@isempty,subpaths(sub_path_row,:)))));

% NB! cellfun for cfg.dataset defines 2:highest populated column
cfg.dataset             = sub_data;
cfg.trialdef.prestim    = 0.500;        % seconds before trigger
cfg.trialdef.poststim   = 0.500;        % seconds after trigger
cfg.trialdef.eventvalue = eval(['temp_stim.' temp_cond{ii} 'trig']) % :/
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

    %offset PO-trials with 50ms in TinMEG1
    %to compensate for ERROR in presentation-script that had trigger wait for gap duration (50ms) even when no gap
    if sub_date.Exp{i} == 'tinmeg1'
    
        for r = 1:length(cfg.trl)
            if (ismember(cfg.trl(r,4), [temp_stim.PO70trig temp_stim.PO60trig]))
                cfg.trl(r,1) = cfg.trl(r,1) - 250; % 250 samples = 50 ms at 5kHz fs
                cfg.trl(r,2) = cfg.trl(r,2) - 250;
                cfg.trl(r,3) = cfg.trl(r,3) - 250;
            end
        end
    end
    
    clear r

% preprocessing
cfg.demean     = 'no';
cfg.lpfilter   = 'no';
cfg.hpfilter   = 'no';
cfg.dftfilter  = 'no';
cfg.channel    = {'MEG', 'ECG', 'EOG'};

preproc = ft_preprocessing(cfg);

    %adjust time-variable to match offset for PO-trials
    if sub_date.Exp{i} == 'tinmeg1'
        if any(ismember(temp_cond{ii}, {'PO60', 'PO70'}))
            cfg = [];
            cfg.offset = 250;
            preproc = ft_redefinetrial(cfg, preproc);   
        end
    end

save([outdir fname], 'preproc', '-v7.3')

%downsample and save
cfg = [];
cfg.resamplefs = fs_ds;

preproc_ds = ft_resampledata(cfg, preproc);

%NB new fname
fname = [temp_cond{ii} '_ds' '.mat'];

save([outdir fname], 'preproc_ds')

%clear temp variables
clear preproc preproc_ds

%end for conditions
end
