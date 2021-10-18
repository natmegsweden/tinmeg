%% Beamformer LCMV
% https://github.com/natmegsweden/meeg_course/blob/master/tutorial_05_beamformer.md

for i = 1:5%length(sub_date.ID);
    
    meg_inpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    mri_inpath = ['../mat_data/MRI_mat/ID' sub_date.ID{i} '/'];
    outdir = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/'];
    
    %load headmodel
    headmodel_meg = load([mri_inpath 'meg_headmodel.mat']); 
    headmodel_meg = headmodel_meg.headmodel_meg;
    
    %Load subject sourcemodel in template grid format (based on MNI)
    subject_grid = load([mri_inpath 'subject_grid.mat']);
    subject_grid = subject_grid.subject_grid;
    
    %Check if subject dir exist, create/define
    if ~exist(outdir, 'file');
    mkdir(outdir);
    end

    %Load all MEG-files
    GOica = load([meg_inpath 'GOica.mat']);
    GOica = GOica.GOica;
    PO60ica = load([meg_inpath 'PO60ica.mat']);
    PO60ica = PO60ica.PO60ica;
    PO70ica = load([meg_inpath 'PO70ica.mat']);
    PO70ica = PO70ica.PO70ica;
    GP60ica = load([meg_inpath 'GP60ica.mat']);
    GP60ica = GP60ica.GP60ica;
    GP70ica = load([meg_inpath 'GP70ica.mat']);
    GP70ica = GP70ica.GP70ica;
    
    %Append data
    cfg = [];
    cfg.keepsampleinfo = 'no'; %if keeping, error because of overlaps
    appended = ft_appenddata(cfg, PO60ica, PO70ica, GP60ica, GP70ica, GOica);
    
    %clear GOica PO60ica PO70ica GP60ica GP70ica;
    
    %Loose ECG/EOG channels
    cfg = [];
    cfg.channel = 'meg';
    appended = ft_selectdata(cfg, appended);
    
    
    %Load or create leadfield
    if ~exist([outdir 'leadfield.mat'], 'file');

    cfg.senstype        = 'meg'; %??
    cfg.grad            = appended.grad;
    cfg.headmodel       = headmodel_meg;
    cfg.sourcemodel     = subject_grid;
    cfg.channel         = 'meg';
    %cfg.grid.resolution = 1;            % Grid spacing 1x1x1 of unit defined below
    %cfg.grid.unit       = 'cm';         % Grid unit

    leadfield = ft_prepare_leadfield(cfg);
    save([outdir 'leadfield.mat'], 'leadfield');
    
    elseif exist([outdir 'leadfield.mat'], 'file');
    
    leadfield = load([outdir 'leadfield.mat']);
    leadfield = leadfield.leadfield;
    
    end
    
    
    %Create noise covarmatrix for denoise_whiten
    cfg.latency = [-0.500 -0.300];
    baseline_noise = ft_selectdata(cfg, appended);
    
    cfg            = [];
    cfg.covariance = 'yes';
    baseline_noise   = ft_timelockanalysis(cfg, baseline_noise);
    
    %Select mags and grads
    selmag  = ft_chantype(baseline_noise.label, 'megmag');
    selgrad = ft_chantype(baseline_noise.label, 'megplanar');
    
    %Denoise_Whiten, see: https://www.fieldtriptoolbox.org/workshop/paris2019/handson_sourceanalysis/
    %Detect the location of the first large 'cliff' in the singular value spectrum of the grads and mags
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
    appended_pw      = ft_denoise_prewhiten(cfg, appended, baseline_noise);
 
    clear baseline_noise appended;
    clear d_mag kappa_mag d_grad kappa_grad u v s_mag s_grad selmag selgrad;
    
    
    %Calculate Filter covariance matrix
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = 'all';
    cfg.channel             = 'MEG';
    data_cov = ft_timelockanalysis(cfg, appended_pw);

    
    %Do initial source analysis to calculte filters
    cfg = [];
    cfg.method              = 'lcmv';
    cfg.channel             = 'meg';
    cfg.lcmv.keepfilter     = 'yes';
    cfg.lcmv.fixedori       = 'yes';
    cfg.lcmv.lambda         = '5%';
    cfg.lcmv.kappa          = kappa;
    cfg.lcmv.projectmom     = 'yes';

    cfg.lcmv.weightnorm     = 'unitnoisegain'; %experiment with this one

    % Original
    cfg.headmodel           = headmodel_meg;
    cfg.sourcemodel         = leadfield;
    source_org = ft_sourceanalysis(cfg, data_cov);
    
    save([outdir 'source_org.mat'], 'source_org');
    
    
    for iii = [1 3]%1:numel(conditions);
    
        %Currently only for PO60, consider: cond.([conditions{1} 'trig'])(1)
        for ii = 1:length(cond.([conditions{iii} 'trig']));

        %Select event data (manual trigger, should refer to structure: cond)
        trigger = cond.([conditions{iii} 'trig'])(ii);

        cfg = [];
        cfg.trials = appended_pw.trialinfo == trigger;
        cfg.latency = [0 0.300];

        stim = ft_selectdata(cfg, appended_pw);

        %Select baseline
        cfg = [];
        cfg.trials = appended_pw.trialinfo == trigger;
        
        if ismember(conditions{iii}, {'GO60', 'GO70', 'PO60', 'PO70'});
        %Baseline window: 200ms before pulse onset in PO trials
        cfg.latency = [-0.200 0];
        elseif ismember(conditions{iii}, {'GP60', 'GP70'});
            %Baseline window variable timepoint: 200ms before gap onset in GP trials
            if ii == 1 %ISI 0
            cfg.latency = [-0.250 -0.050];
            elseif ii == 2 %ISI 60
                cfg.latency = [-0.310 -0.110];
            elseif ii == 3 %ISI 120
                cfg.latency = [-0.370 -0.170];
            elseif ii == 4 %ISI 240
                cfg.latency = [-0.490 -0.290];
            end
        end

        base = ft_selectdata(cfg, appended_pw);


        %Compute covariance matrix for data and baseline
        cfg = [];
        cfg.covariance              = 'yes';
        cfg.covariancewindow        = 'all';
        cfg.preproc.demean          = 'yes';
        cfg.keeptrials              = 'no'; %Y/N makes no difference?

        stim_cov = ft_timelockanalysis(cfg, stim);
        base_cov = ft_timelockanalysis(cfg, base);
        
        clear stim base;

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

        save([outdir cond.([conditions{iii} 'label']){ii} '_stim_source.mat'], 'stim_source');
        save([outdir cond.([conditions{iii} 'label']){ii} '_base_source.mat'], 'base_source');
        
        cfg = [];
        cfg.parameter = 'pow';
        cfg.operation = 'x1 - mean(x2)'; %Different window duration
        powdiff_sub = ft_math(cfg, stim_source, base_source);
        
        save([outdir cond.([conditions{iii} 'label']){ii} '_powdiff.mat'], 'powdiff_sub');
    
        %For trial
        end
        
    %For condition
    end

%For subject
end
