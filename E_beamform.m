%% https://github.com/natmegsweden/meeg_course/blob/master/tutorial_05_beamformer.md


%% Load template grid

% Loads in "Full_analysis.m" for consistency with MR preprocessing.
% load('/../../fieldtrip-20210311/template/sourcemodel/standard_sourcemodel3d10mm');
% template_grid = sourcemodel;
% clear sourcemodel

%% inspect subjects headmodel fit to grid

for i = 1%:4%length(sub_date.ID);
    
    inpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    outdir = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/'];
    
    %headmodel_meg
    headmodel_meg = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '_MEG_headmodel.mat']); 
    headmodel_meg = headmodel_meg.headmodel_meg;
    
    subject_grid = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '_sub_grid.mat']);
    subject_grid = subject_grid.subject_grid;
    
    
%     %Final plot - aligned MEG   
%     fig = figure;
%     hold on;
%     ft_plot_headmodel(headmodel_meg, 'edgecolor', 'none', 'facealpha', 0.4,'facecolor', 'b');
%     ft_plot_mesh(subject_grid.pos(subject_grid.inside,:));
%     title(['SUBJECT: ' sub_date.ID{i}]);
% 
%     %Pause loop until figure is closed
%     uiwait(fig);
    

    %Check if subject dir exist, create/define
    if ~exist(outdir, 'file');
    mkdir(outdir);
    end

    %Load in the stupidest way possible
    GOica = load([inpath 'GOica.mat']);
    GOica = GOica.GOica;
    PO60ica = load([inpath 'PO60ica.mat']);
    PO60ica = PO60ica.PO60ica;
    PO70ica = load([inpath 'PO70ica.mat']);
    PO70ica = PO70ica.PO70ica;
    GP60ica = load([inpath 'GP60ica.mat']);
    GP60ica = GP60ica.GP60ica;
    GP70ica = load([inpath 'GP70ica.mat']);
    GP70ica = GP70ica.GP70ica;

    
    %Append data
    cfg = [];
    cfg.keepsampleinfo = 'no'; %if keeping, error because of overlaps
    appended = ft_appenddata(cfg, PO60ica, PO70ica, GP60ica, GP70ica, GOica);
    
    
    cfg.senstype        = 'meg'; %??
    cfg.grad            = appended.grad;
    cfg.headmodel       = headmodel_meg;
    cfg.sourcemodel     = subject_grid;
    cfg.channel         = 'meg';
    % cfg.grid.resolution = 1;            % Grid spacing 1x1x1 of unit defined below
    % cfg.grid.unit       = 'cm';         % Grid unit

    leadfield_meg = ft_prepare_leadfield(cfg);
    
    save([outdir 'leadfield.mat'], 'leadfield_meg');
    
    
    
    
    %Select filter data (manual trigger, should refer to structure: cond)
    cfg = [];
    cfg.trials = appended.trialinfo == 33032;
    cfg.latency = [-300.0 0.300];

    PO60_95_filter = ft_selectdata(cfg, appended);
        
    %Select event data (manual trigger, should refer to structure: cond)
    cfg = [];
    cfg.trials = appended.trialinfo == 33032;
    cfg.latency = [0.0 0.300];

    PO60_95 = ft_selectdata(cfg, appended);
    
    %Select baseline
    cfg = [];
    cfg.trials = appended.trialinfo == 33032;
    cfg.latency = [-300 0];
    
    PO60_95_base = ft_selectdata(cfg, appended);
    
    
    
    %Calculate Filter covariance matrix and Kappa (rank deficiency)
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = 'all';
    cfg.channel             = 'MEG';
    data_cov = ft_timelockanalysis(cfg, PO60_95_filter);

    [u,s,v] = svd(data_cov.cov);
    d       = -diff(log10(diag(s)));
    d       = d./std(d);
    kappa   = find(d>5,1,'first');
    fprintf('Kappa = %i\n', kappa)

    %consider rank(data_cov.cov)?

    % figure;
    % semilogy(diag(s),'o-');   
    
    
    %Do initial source analysis to calculte filters
    cfg = [];
    cfg.method              = 'lcmv';
    cfg.channel             = 'meg';
    cfg.lcmv.keepfilter     = 'yes';
    cfg.lcmv.fixedori       = 'yes';
    cfg.lcmv.lambda         = '5%';
    cfg.lcmv.kappa          = kappa;
    cfg.lcmv.projectmom     = 'yes';

    % Original
    cfg.headmodel           = headmodel_meg;
    cfg.sourcemodel         = leadfield_meg;
    source_org = ft_sourceanalysis(cfg, data_cov);
    
    save([outdir 'source_org.mat'], 'source_org');
    


    %Compute covariance matrix for data and baseline
    cfg = [];
    cfg.covariance              = 'yes';
    cfg.covariancewindow        = 'all';
    cfg.preproc.demean          = 'yes';
    cfg.keeptrials              = 'yes'; %as R/L AC is highly correlated

    PO60_95_cov = ft_timelockanalysis(cfg, PO60_95);
    PO60_95_base_cov = ft_timelockanalysis(cfg, PO60_95_base);
    
    
    %Source analysis on baseline and stim data
    cfg=[];
    cfg.method              = 'lcmv';
    cfg.grid                = leadfield_meg;
    cfg.sourcemodel.filter  = source_org.avg.filter;  % Reuse avg filter
    cfg.headmodel           = headmodel_meg;
    cfg.channel             = 'meggrad';
    cfg.senstype            = 'MEG';

    PO60_95_base_source = ft_sourceanalysis(cfg, PO60_95_base_cov);
    PO60_95_source = ft_sourceanalysis(cfg, PO60_95_cov);
    
    %Contrast between stim and baseline
    contrast_lcmv = PO60_95_source;       % Copy
    contrast_lcmv.avg.pow = (PO60_95_source.avg.pow-PO60_95_base_source.avg.pow)./PO60_95_base_source.avg.pow;
    
    mri_resliced = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '_mri_resliced.mat']);
    mri_resliced = mri_resliced.mri_resliced;

    %Interpolate and plot contrasted source
    cfg = [];
    cfg.parameter    = 'pow';
    cfg.interpmethod = 'nearest';
    source_int  = ft_sourceinterpolate(cfg, contrast_lcmv, mri_resliced);
    
    save([outdir 'PO6095_interpolated.mat'], 'source_int');
    
end

%Plot
cfg = [];
cfg.method          = 'ortho';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';
ft_sourceplot(cfg, source_int);

