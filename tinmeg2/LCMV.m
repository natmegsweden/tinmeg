%% Source reconstruction (LCMV Beamformer)

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

%Create Leadfield
cfg = [];
cfg.senstype            = 'MEG';
cfg.grad                = data.grad;
cfg.headmodel           = headmodel;
cfg.channel             = 'MEG';
cfg.grid.resolution     = 1;
cfg.grid.unit           = 'cm';

leadfield = ft_prepare_leadfield(cfg);

%Create timelockeds for trial of interest
cfg = [];
cfg.trials = data.trialinfo == 4; % 4 is PO trigger
cfg.covariance = 'yes';
cfg.covariancewindow = 'all';
cfg.preproc.demean = 'yes';
full_trial = ft_timelockanalysis(cfg, data);

%
cfg = [];
cfg.latency = [-0.500 -0.250];

baseline_epoch = ft_selectdata(cfg, full_trial);

cfg.latency = [0 0.500];
toi_epoch =  ft_selectdata(cfg, full_trial);

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

cfg.sourcemodel.filter = source_all.avg.filter

source_baseline = ft_sourceanalysis(cfg, baseline_epoch);
source_toi = ft_sourceanalysis(cfg, toi_epoch);

%Make contrast
lcmv_contrast = source_all; %copy
lcmv_contrast.avg.pow = (source_toi.avg.pow - source_baseline.avg.pow) ./ source_baseline.avg.pow;

%Interpolate contrast on MRI
cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';

source_int = ft_sourceinterpolate(cfg, lcmv_contrast, mri_resliced)

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolormap = 'jet';

ft_sourceplot(cfg, source_int);



