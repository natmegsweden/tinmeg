%% Source reconstruction (LCMV Beamformer)

%Try ft_timelockeds separately for baseline_epoch & toi_epoch

addpath('../');

%Load MRI for interpolation/plot
mri_resliced = load('../mat_data/MRI_mat_tinmeg2/ID0916/mri_resliced.mat');
mri_resliced = mri_resliced.mri_resliced;

mri_resliced = ft_convert_units(mri_resliced, 'cm');

%Load timelocked data
data = load('../mat_data/ICA/ID0916/tin0_bkg0_ica.mat');
data = data.tin0_bkg0_ica;

%load headmodel
headmodel = load('../mat_data/MRI_mat_tinmeg2/ID0916/meg_headmodel.mat');
headmodel = headmodel.headmodel_meg;

headmodel = ft_convert_units(headmodel, 'cm');

%Create Leadfield
cfg = [];
cfg.senstype            = 'MEG';
cfg.grad                = data.grad;
cfg.headmodel           = headmodel;
cfg.channel             = 'MEG';
cfg.resolution          = 1;

leadfield = ft_prepare_leadfield(cfg);

%Select data
cfg = [];
cfg.trials = data.trialinfo == 4; % 4 is PO trigger
cfg.latency = [-0.500 0.500];

full_epoch = ft_selectdata(cfg, data);

cfg.latency = [-0.500 0];
baseline_epoch = ft_selectdata(cfg, data);

cfg.latency = [0 0.500];
toi_epoch = ft_selectdata(cfg, data);

% Timelockeds for full epoch
cfg = [];
cfg.covariance         = 'yes';
cfg.covariancewindow   = 'all';
cfg.preproc.demean     = 'yes';
full_trial = ft_timelockanalysis(cfg, full_epoch);

% Timelockeds fo baseline and stimulation
cfg = [];
cfg.covariance         = 'yes';
cfg.covariancewindow   = 'all';
baseline_epoch = ft_timelockanalysis(cfg, baseline_epoch);
toi_epoch = ft_timelockanalysis(cfg, toi_epoch);

%Source reconstruction on full_trial to construct common spatial filter
cfg = [];
cfg.method = 'lcmv';
cfg.grid = leadfield;
cfg.headmodel = headmodel;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.lambda = '5%';
cfg.channel = 'MEG';
cfg.senstype = 'MEG';

source_all = ft_sourceanalysis(cfg, full_trial);

cfg.sourcemodel.filter = source_all.avg.filter;
cfg.sourcemodel.filterdimord = source_all.avg.filterdimord;

source_baseline = ft_sourceanalysis(cfg, baseline_epoch);
source_toi = ft_sourceanalysis(cfg, toi_epoch);

%Make contrast
lcmv_contrast = source_toi; %copy

lcmv_contrast.avg.pow = (source_toi.avg.pow - source_baseline.avg.pow) ./ source_baseline.avg.pow;

%Interpolate contrast on MRI
cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';

source_int = ft_sourceinterpolate(cfg, lcmv_contrast, mri_resliced);

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolormap = 'jet';

ft_sourceplot(cfg, source_int);

