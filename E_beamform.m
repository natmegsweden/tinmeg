
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

load(['../mat_data/MRI_mat/ID' sub_date.ID{2} '_MEG_headmodel.mat'])

load(['../mat_data/ID' sub_date.ID{2} '_PO60_ds_clean.mat'])

%put loaded data in structure?


%% Make leadfields for MEG: magnetometers
cfg.senstype        = 'meg';
cfg.grad            = cleaned4mat.grad;
cfg.headmodel       = headmodel_meg;
cfg.channel         = 'meg';
cfg.grid.resolution = 1;            % Grid spacing 1x1x1 of unit defined below
cfg.grid.unit       = 'cm';         % Grid unit

leadfield_meg = ft_prepare_leadfield(cfg);

%% Select data
cfg = [];
cfg.trials = cleaned4mat.trialinfo == 33032;
cfg.latency = [0.100 0.200];

data = ft_selectdata(cfg, cleaned4mat);

%% evoked
cfg = [];
cfg.covariance              = 'yes';
cfg.covariancewindow        = 'all';
cfg.preproc.demean          = 'yes';

evoked = ft_timelockanalysis(cfg, data);

%% LCMV
cfg = [];
cfg.method          = 'lcmv';
cfg.grid            = leadfield_meg;
cfg.headmodel       = headmodel_meg;
cfg.lcmv.lambda     = '5%';
cfg.channel         = 'meggrad';
cfg.senstype        = 'MEG';
source_lcmv = ft_sourceanalysis(cfg, evoked);

%% Load MRI
load(['../mat_data/MRI_mat/ID' sub_date.ID{2} '_mri_resliced.mat']);

%% Interpolate
cfg = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, source_lcmv, mri_resliced);

%% Plot
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
ft_sourceplot(cfg,source_int);

