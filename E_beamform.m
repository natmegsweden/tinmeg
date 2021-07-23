%% https://github.com/natmegsweden/meeg_course/blob/master/tutorial_05_beamformer.md


%% Load template grid

% Loads in "Full_analysis.m" for consistency with MR preprocessing.
% load('/../../fieldtrip-20210311/template/sourcemodel/standard_sourcemodel3d10mm');
% template_grid = sourcemodel;
% clear sourcemodel

%% inspect subjects headmodel fit to grid

for i = 1:4%length(sub_date.ID);
    
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

%     %Load in the stupidest way possible
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
    
    
    %Load or create leadfield
    if ~exist([outdir 'leadfield.mat'], 'file');

    cfg.senstype        = 'meg'; %??
    cfg.grad            = appended.grad;
    cfg.headmodel       = headmodel_meg;
    cfg.sourcemodel     = subject_grid;
    cfg.channel         = 'meg';
    % cfg.grid.resolution = 1;            % Grid spacing 1x1x1 of unit defined below
    % cfg.grid.unit       = 'cm';         % Grid unit

    leadfield = ft_prepare_leadfield(cfg);
    save([outdir 'leadfield.mat'], 'leadfield');
    
    elseif exist([outdir 'leadfield.mat'], 'file');
    
    leadfield = load([outdir 'leadfield.mat']);
    leadfield = leadfield.leadfield;
    
    end
    
    %Loose ECG/EOG channels
%     cfg = [];
%     cfg.channel = 'meg';
%     PO60ica = ft_selectdata(cfg, PO60ica);
    
    %Create noise covarmatrix for denoise_whiten
    cfg.latency = [-0.350 -0.100];
    baseline_noise = ft_selectdata(cfg, appended);
    
    cfg            = [];
    cfg.covariance = 'yes';
    baseline_noise   = ft_timelockanalysis(cfg, baseline_noise); 
    
    %to selects mags and grads
    selmag  = ft_chantype(baseline_noise.label, 'megmag');
    selgrad = ft_chantype(baseline_noise.label, 'megplanar');
    
    %Denoise_Whiten
    % the following lines detect the location of the first large 'cliff' in the singular value spectrum of the grads and mags
    [u,s_mag,v]  = svd(baseline_noise.cov(selmag,  selmag));
    [u,s_grad,v] = svd(baseline_noise.cov(selgrad, selgrad));
    d_mag = -diff(log10(diag(s_mag))); d_mag = d_mag./std(d_mag);
    kappa_mag = find(d_mag>4,1,'first');
    d_grad = -diff(log10(diag(s_grad))); d_grad = d_grad./std(d_grad);
    kappa_grad = find(d_grad>4,1,'first');
    
    kappa = min(kappa_mag,kappa_grad);
    
    cfg            = [];
    cfg.channel    = 'meg';
    cfg.kappa      = kappa;
    dataw_meg      = ft_denoise_prewhiten(cfg, appended, baseline_noise);
 
%     cfg.layout = 'neuromag306all.lay';
%     ft_multiplotER(cfg, PO60ica);
    
    
    %Calculate Filter covariance matrix and Kappa (rank deficiency)
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = 'all';
    cfg.channel             = 'MEG';
    data_cov = ft_timelockanalysis(cfg, dataw_meg);

%     [u,s,v] = svd(data_cov.cov);
%     d       = -diff(log10(diag(s)));
%     d       = d./std(d);
%     kappa   = find(d>5,1,'first');
%     fprintf('Kappa = %i\n', kappa)

    %consider rank(data_cov.cov)?

    % figure;
    % semilogy(diag(s),'o-');  
    
    %ft_inverse_lcmv
    
    %Do initial source analysis to calculte filters
    cfg = [];
    cfg.method              = 'lcmv';
    cfg.channel             = 'meg';
    cfg.keepfilter          = 'yes';
    cfg.lcmv.fixedori       = 'yes';
    cfg.lambda              = '5%';
    cfg.kappa               = kappa;
    cfg.lcmv.projectmom     = 'yes';

    cfg.lcmv.weightnorm     = 'unitnoisegain'; %experiment with this one

    % Original
    cfg.headmodel           = headmodel_meg;
    cfg.sourcemodel         = leadfield;
    source_org = ft_sourceanalysis(cfg, data_cov);
    
    save([outdir 'source_org.mat'], 'source_org');
    
    mri_resliced = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '_mri_resliced.mat']);
    mri_resliced = mri_resliced.mri_resliced;
    
    for ii = 1:length(cond.PO60trig);
    %Select filter data (manual trigger, should refer to structure: cond)
    trigger = cond.PO60trig(ii);
        
    
    %Select event data (manual trigger, should refer to structure: cond)
    cfg = [];
    cfg.trials = dataw_meg.trialinfo == trigger;
    cfg.latency = [0 0.300];

    stim = ft_selectdata(cfg, dataw_meg);
    
    %Select baseline
    cfg = [];
    cfg.trials = dataw_meg.trialinfo == trigger;
    cfg.latency = [-0.300 0];
    
    base = ft_selectdata(cfg, dataw_meg);
    

    %Compute covariance matrix for data and baseline
    cfg = [];
    cfg.covariance              = 'yes';
    cfg.covariancewindow        = 'all';
    cfg.preproc.demean          = 'yes';
    cfg.keeptrials              = 'no'; %Y/N makes no difference

    stim_cov = ft_timelockanalysis(cfg, stim);
    base_cov = ft_timelockanalysis(cfg, base);
    
    
    %Source analysis on baseline and stim data
    cfg=[];
    cfg.method              = 'lcmv';
    cfg.sourcemodel         = leadfield;
    cfg.sourcemodel.filter  = source_org.avg.filter;  % Reuse avg filter
    cfg.headmodel           = headmodel_meg;
    cfg.channel             = 'meg'; %grad & mag
    cfg.senstype            = 'MEG';

    stim_source = ft_sourceanalysis(cfg, stim_cov);
    base_source = ft_sourceanalysis(cfg, base_cov);
    
    
    %Contrast between stim and baseline
    contrast_lcmv = stim_source;       % Copy
    contrast_lcmv.avg.pow = (stim_source.avg.pow-base_source.avg.pow)./base_source.avg.pow;


    %Interpolate contrasted source
    cfg = [];
    cfg.parameter    = 'pow';
    cfg.interpmethod = 'nearest';
    source_int  = ft_sourceinterpolate(cfg, contrast_lcmv, mri_resliced);
    
    save([outdir char(cond.PO60label(ii)) '_interpolated.mat'], 'source_int');
    
    %Plot and save
    cfg = [];
    cfg.method          = 'ortho';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';
    
%     cfg.maskparameter = cfg.funparameter
%     cfg.colorlim      = [0 3] % or 'zeromax'
%     cfg.opacitymap    = 'rampup'
%     cfg.opacitylim    = [0 3] % or 'zeromax'

    
    cfg.position        = [700 300 950 950];
    ft_sourceplot(cfg, source_int);

    saveas(gcf, [outdir char(cond.PO60label(ii)) '.png']);
    
    close
    
    end

    
end

