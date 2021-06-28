%% https://github.com/natmegsweden/meeg_course/blob/master/tutorial_05_beamformer.md


%% Load template grid

load('/../../fieldtrip-20210311/template/sourcemodel/standard_sourcemodel3d10mm');
template_grid = sourcemodel;
clear sourcemodel

%% Load subject data
load(['../mat_data/MRI_mat/ID' sub_date.ID{2} '_MEG_headmodel.mat']);  
load(['../mat_data/MRI_mat/ID' sub_date.ID{2} '_mri_resliced.mat']);

%load(['../mat_data/ID' sub_date.ID{2} '_PO60_ds_clean.mat']);

%put loaded data in structure?
%ladda in headmodels + mr_resliced: skapa normaliserad source_model
%spara sourcemodels

%% Create the subject specific grid

cfg           = [];
cfg.method    = 'basedonmni';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri_resliced;
cfg.unit      = 'mm';

subject_grid = ft_prepare_sourcemodel(cfg);

save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) 'sub_grid'], 'subject_grid');

% make a figure of the single subject headmodel, and grid positions
% figure; hold on;
% ft_plot_headmodel(headmodel_meg, 'edgecolor', 'none', 'facealpha', 0.4);
% ft_plot_mesh(grid.pos(grid.inside,:));


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
cfg.latency = [0.0 0.150];

data = ft_selectdata(cfg, cleaned4mat);

%% evoked
cfg = [];
cfg.covariance              = 'yes';
cfg.covariancewindow        = 'all';
cfg.preproc.demean          = 'yes';
cfg.keeptrials              = 'yes'; %as R/L AC is highly correlated

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

%% Interpolate (middle of head bias (?))
cfg = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, source_lcmv, mri_resliced);

% Plot
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
ft_sourceplot(cfg,source_int);


%% Select full window, baseline and TOI

%full window
cfg = [];
cfg.trials = cleaned4mat.trialinfo == 33032;
cfg.latency = [-0.200 0.200];

data_all = ft_selectdata(cfg, cleaned4mat);

%baseline window
cfg.latency = [-0.200 -0.050];

data_base = ft_selectdata(cfg, cleaned4mat);

cfg.latency = [0.0 0.150];

data_stim = ft_selectdata(cfg, cleaned4mat);

%% Evoked for TOI
% Combined
cfg = [];
cfg.covariance         = 'yes';
cfg.covariancewindow   = 'all';
cfg.preproc.demean     = 'yes';
evo_all = ft_timelockanalysis(cfg, data_all);

% Baseline and stimulation
cfg = [];
cfg.covariance         = 'yes';
cfg.covariancewindow   = 'all';
evo_base = ft_timelockanalysis(cfg, data_base);
evo_stim = ft_timelockanalysis(cfg, data_stim);

%% Calculate beamformer filter - average for TOI

cfg=[];
cfg.method          = 'lcmv';
cfg.grid            = leadfield_meg;
cfg.headmodel       = headmodel_meg;
cfg.lcmv.keepfilter = 'yes';        % save the filter in the output data
cfg.lcmv.lambda     = '5%';
cfg.channel         = 'meggrad';
cfg.senstype        = 'MEG';

source_all = ft_sourceanalysis(cfg, evo_all);

%% Source analysis on baseline and stim data
cfg=[];
cfg.method              = 'lcmv';
cfg.grid                = leadfield_meg;
cfg.sourcemodel.filter  = source_all.avg.filter;  % Reuse avg filter
cfg.headmodel           = headmodel_meg;
cfg.channel             = 'meggrad';
cfg.senstype            = 'MEG';

source_base = ft_sourceanalysis(cfg, evo_base);
source_stim = ft_sourceanalysis(cfg, evo_stim);

%% Create contrast - baseline and stim

contrast_lcmv = source_stim;       % Copy
contrast_lcmv.avg.pow = (source_stim.avg.pow-source_base.avg.pow)./source_base.avg.pow;

%Save here

%% Interpolate and plot contrasted source
cfg = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, contrast_lcmv, mri_resliced);

%Plot
cfg = [];
cfg.method          = 'ortho';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';
ft_sourceplot(cfg, source_int);

