%% Beamformer LCMV
% https://github.com/natmegsweden/meeg_course/blob/master/tutorial_05_beamformer.md

%Contrasts of interest
coi = {'POvsB', 'GPPvsB', 'GPPvsPO'};

for i = 10 %1:25%length(sub_date.ID);
    
    meg_inpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    mri_inpath = ['../mat_data/MRI_mat/ID' sub_date.ID{i} '/'];
    outdir = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/'];
    
    %load headmodel
    headmodel_meg = load([mri_inpath 'meg_headmodel.mat']); 
    headmodel_meg = headmodel_meg.headmodel_meg;
    headmodel_meg = ft_convert_units(headmodel_meg, 'cm');
    
    %Load subject sourcemodel in template grid format (based on MNI)
    %subject_grid = load([mri_inpath 'subject_grid.mat']);
    %subject_grid = subject_grid.subject_grid;
    
    %Check if subject dir exist, create/define
    if ~exist(outdir, 'file');
    mkdir(outdir);
    end

    %Load MEG files of interest
    %GOica = load([meg_inpath 'GOica.mat']);
    %GOica = GOica.GOica;
    PO60ica = load([meg_inpath 'PO60ica.mat']);
    PO60ica = PO60ica.PO60ica;
    %PO70ica = load([meg_inpath 'PO70ica.mat']);
    %PO70ica = PO70ica.PO70ica;
    GP60ica = load([meg_inpath 'GP60ica.mat']);
    GP60ica = GP60ica.GP60ica;
    %GP70ica = load([meg_inpath 'GP70ica.mat']);
    %GP70ica = GP70ica.GP70ica;
    
    %Append data
    cfg = [];
    cfg.keepsampleinfo = 'no'; %if keeping, error because of overlaps
    %appended = ft_appenddata(cfg, PO60ica, PO70ica, GP60ica, GP70ica, GOica);
    %appended = ft_appenddata(cfg, PO60ica, GP60ica);
    
    %Specify triggers
    POtrig = 33288;
    GPPtrig = 49688;

    % TEST
    cfg = [];
    cfg.channel = 'meg';
    cfg.trials = PO60ica.trialinfo == POtrig;

    appended = ft_selectdata(cfg, PO60ica);
    
    %Select prestim noise for prewhiten
    cfg =[];
    cfg.latency = [-0.500 -0.250];
    baseline_noise = ft_selectdata(cfg, appended);
    
    %Keep trials of interest and MEG channels only in baseline_noise
    cfg = [];
    cfg.channel = 'meg';
    %PO60_90 = 33288, GP60_i240 = 49688
    cfg.trials = baseline_noise.trialinfo == POtrig | baseline_noise.trialinfo == GPPtrig;
    baseline_noise = ft_selectdata(cfg, baseline_noise);

    %Keep trials of interest and MEG channels only in data
    cfg = [];
    cfg.channel = 'meg';
    %PO60_90 = 33288, GP60_i240 = 49688
    cfg.trials = appended.trialinfo == POtrig | appended.trialinfo == GPPtrig;
    data = ft_selectdata(cfg, appended);

    %clear GOica PO60ica PO70ica GP60ica GP70ica;
   
    %Create covmatrix for denoise_whiten using baseline_noise
    cfg            = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'all';
    baseline_noise   = ft_timelockanalysis(cfg, baseline_noise);
    
    %Select mags and grads
    selmag  = ft_chantype(baseline_noise.label, 'megmag');
    selgrad = ft_chantype(baseline_noise.label, 'megplanar');
    
    % C = baseline_noise.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
    % figure;imagesc(C);hold on;plot(102.5.*[1 1],[0 306],'w','linewidth',2);plot([0 306],102.5.*[1 1],'w','linewidth',2);

    %Denoise_Whiten, see: https://www.fieldtriptoolbox.org/workshop/paris2019/handson_sourceanalysis/
    %Detect the location of the first large 'cliff' in the singular value spectrum of the grads and mags
    [u,s_mag,v]  = svd(baseline_noise.cov(selmag,  selmag));
    [u,s_grad,v] = svd(baseline_noise.cov(selgrad, selgrad));
    d_mag = -diff(log10(diag(s_mag))); d_mag = d_mag./std(d_mag);
    kappa_mag = find(d_mag>4,1,'first');
    d_grad = -diff(log10(diag(s_grad))); d_grad = d_grad./std(d_grad);
    kappa_grad = find(d_grad>4,1,'first');
    
    kappa = min(kappa_mag,kappa_grad);
    
    % Singular values of mag and grad covariance matrix
    % figure;plot(log10(diag(s_mag)),'o');
    % figure;plot(log10(diag(s_grad)),'o');
    
    cfg            = [];
    cfg.channel    = 'meg';
    cfg.kappa      = kappa;
    data_pw    = ft_denoise_prewhiten(cfg, data, baseline_noise);

    save([outdir 'pw_tlks.mat'], "data_pw");

    %clear baseline_noise appended;
    clear d_mag kappa_mag d_grad kappa_grad u v s_mag s_grad selmag selgrad;
    
    %Load or create leadfield
    if ~exist([outdir 'leadfield.mat'], 'file');

    cfg.senstype        = 'meg';
    cfg.grad            = data_pw.grad;
    cfg.headmodel       = headmodel_meg;
    %cfg.sourcemodel     = subject_grid;
    cfg.channel         = 'meg';
    cfg.grid.resolution = 1;
    cfg.grid.unit       = 'cm';

    leadfield = ft_prepare_leadfield(cfg);
    save([outdir 'leadfield.mat'], 'leadfield');
    
    elseif exist([outdir 'leadfield.mat'], 'file');
    
    leadfield = load([outdir 'leadfield.mat']);
    leadfield = leadfield.leadfield;
    
    end

    %Covariance for full epoch all conditions
    cfg = [];
    cfg.preproc.demean = 'yes';
    %cfg.preproc.baselinewindow = [-0.500 -0.250];
    cfg.covariance = 'yes';
    data_cov = ft_timelockanalysis(cfg, data_pw);

%     C = data_cov.cov([find(selmag);find(selgrad)],[find(selmag);find(selgrad)]);
%     figure;imagesc(C);hold on;plot(102.5.*[1 1],[0 306],'w','linewidth',2);plot([0 306],102.5.*[1 1],'w','linewidth',2);

    %Source reconstruction on full_trial to construct common spatial filter
    cfg = [];
    cfg.method = 'lcmv';
    cfg.grid = leadfield;
    cfg.headmodel = headmodel_meg;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.lambda = '5%';
    cfg.channel = 'MEG';
    cfg.senstype = 'MEG';
    
    source_all = ft_sourceanalysis(cfg, data_cov);
    
    %For contrasts of interest
    for ii = 1:numel(coi)

        %Testing the basic one: PulseOnly response vs pre-stim baseline
        if ismember(coi(ii), 'POvsB')

            %Select baseline and toi
            cfg = [];
            cfg.trials = data_pw.trialinfo == POtrig;
            cfg.latency = [-0.500 -0.250]; % -200 to 0 ms before pulse
            baseline = ft_selectdata(cfg, data_pw);

            %Select experimental stim condition and toi
            cfg = [];
            cfg.trials = data_pw.trialinfo == POtrig;
            cfg.latency = [0.050 0.150]; % 50 to 150 ms after pulse
            stim = ft_selectdata(cfg, data_pw);
        
            % Timelockeds for baseline and stimulation
            cfg = [];
            cfg.covariance         = 'yes';
            cfg.covariancewindow   = 'all';
            baseline_cov = ft_timelockanalysis(cfg, baseline);
            stim_cov = ft_timelockanalysis(cfg, stim);
            
            %config for sourcemodel stim vs baseline
            cfg = [];
            cfg.method = 'lcmv';
            cfg.grid = leadfield;
            cfg.headmodel = headmodel_meg;
            cfg.lcmv.keepfilter = 'yes';
            cfg.lcmv.lambda = '5%';
            cfg.channel = 'MEG';
            cfg.senstype = 'MEG';
            %Use sourcemodel.filter for all full trials
            cfg.sourcemodel.filter = source_all.avg.filter;
            cfg.sourcemodel.filterdimord = source_all.avg.filterdimord;
            
            source_baseline = ft_sourceanalysis(cfg, baseline_cov);
            source_stim = ft_sourceanalysis(cfg, stim_cov);
            
            %Make contrast
            lcmv_contrast = source_stim; %copy

            %invert GPPvsPO so "inhibition" is positive
            if ismember(coi(ii), 'GPPvsPO'); 
                lcmv_contrast.avg.pow = (source_baseline.avg.pow - source_stim.avg.pow) ./ source_baseline.avg.pow;
            else
                lcmv_contrast.avg.pow = (source_stim.avg.pow - source_baseline.avg.pow) ./ source_baseline.avg.pow;
            end

            %Interpolate contrast on template MRI
            cfg = [];
            cfg.parameter = 'pow';
            cfg.interpmethod = 'nearest';
            
            load([mri_inpath 'mri_resliced.mat'])

            source_int = ft_sourceinterpolate(cfg, lcmv_contrast, mri_resliced);

            figure('Position', [683,74,1095,1091])
            cfg = [];
            cfg.method = 'ortho';
            cfg.funparameter = 'pow';
            cfg.funcolormap = 'jet';

            ft_sourceplot(cfg, source_int);

            %saveas(gcf, ['../mat_data/source_reconstruction/testfigures2/' sub_date.ID{i} '.png'])
            %close;

        end
    
    %for contrasts
    end

%For subject
end
