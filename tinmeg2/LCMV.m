%% Source reconstruction (LCMV Beamformer) of PO and GPP

%Check agreement of PO/GP-plots per sensor and 'inhibition vector'
% 'MinProminence' in islocalmin

%load tois from sensor analysis
load('../mat_data/timelockeds/tinmeg2/tois.mat');

addpath('../');

% disp([time(tois.tin0_bkg0_1{1,1}) time(tois.tin0_bkg0_1{1,2})]);
% disp([time(tois.tin0_bkg3_1{1,1}) time(tois.tin0_bkg3_1{1,2})]);
% disp([time(tois.tin0_bkg8_1{1,1}) time(tois.tin0_bkg8_1{1,2})]);
% 
% disp([time(tois.tin0_bkg0_2{1,1}) time(tois.tin0_bkg0_2{1,2})]);
% disp([time(tois.tin0_bkg3_2{1,1}) time(tois.tin0_bkg3_2{1,2})]);
% disp([time(tois.tin0_bkg8_2{1,1}) time(tois.tin0_bkg8_2{1,2})]);

bkgs = [0 3 8];

plots = {'POvsB', 'GPPvsB', 'GPPvsPO'};

%For subject
for i = 1:numel(sub_date.ID)

        %Load MRI for interpolation/plot
        mri_resliced = load(['../mat_data/MRI_mat_tinmeg2/ID' sub_date.ID{i} '/mri_resliced.mat']);
        mri_resliced = mri_resliced.mri_resliced;
        
        mri_resliced = ft_convert_units(mri_resliced, 'cm');

        %load headmodel
        headmodel = load(['../mat_data/MRI_mat_tinmeg2/ID' sub_date.ID{i} '/meg_headmodel.mat']);
        headmodel = headmodel.headmodel_meg;
        
        headmodel = ft_convert_units(headmodel, 'cm');
    
    for ii = 1:numel(bkgs)
            
            %condition name
            name = ['tin0_bkg' num2str(bkgs(ii))];
    
            %Load timelocked data
            data = load(['../mat_data/ICA/ID' sub_date.ID{i} '/' name '_ica.mat']);
            data = data.([name '_ica']);
            
            %Get time vector from data
            time = data.time{1,1};
        
            %Create Leadfield
            cfg = [];
            cfg.senstype            = 'MEG';
            cfg.grad                = data.grad;
            cfg.headmodel           = headmodel;
            cfg.channel             = 'MEG';
            cfg.resolution          = 1;
            
            leadfield = ft_prepare_leadfield(cfg);

            POtrig = cond.([name 'trig'])(3);
            GPPtrig = cond.([name 'trig'])(1);
            
            for iii = 1:numel(plots);
                %Select data
                cfg = [];
                
                if ismember(plots(iii), 'POvsB'); 
                    %Select PO-trials
                    cfg.trials = data.trialinfo == POtrig;
                    cfg.latency = [-0.500 0.500];
                    full_epoch = ft_selectdata(cfg, data);

                    %Select baseline
                    cfg.latency = [-0.500 -0.250];
                    baseline_epoch = ft_selectdata(cfg, data);
                    
                    %Select toi
                    cfg.latency = [0.050 0.150];
                    toi_epoch = ft_selectdata(cfg, data);
                
                elseif ismember(plots(iii), 'GPPvsB'); cfg.trials = data.trialinfo == GPPtrig;
                    %Select GPP-trials
                    cfg.trials = data.trialinfo == GPPtrig;
                    cfg.latency = [-0.500 0.500];
                    full_epoch = ft_selectdata(cfg, data);

                    %Select baseline
                    cfg.latency = [-0.500 -0.250];
                    baseline_epoch = ft_selectdata(cfg, data);
                    
                    %Select toi
                    cfg.latency = [0.050 0.150];
                    toi_epoch = ft_selectdata(cfg, data);
                
                elseif ismember(plots(iii), 'GPPvsPO'); 
                    %Select GPP AND PO trials
                    cfg.trials = data.trialinfo == POtrig | data.trialinfo == GPPtrig;
                    cfg.latency = [-0.500 0.500];
                    full_epoch = ft_selectdata(cfg, data);

                    %Select baseline
                    cfg.trials = data.trialinfo == POtrig %PO toi is baseline
                    cfg.latency = [0.050 0.150];
                    baseline_epoch = ft_selectdata(cfg, data);
                    
                    cfg.trials = data.trialinfo == GPPtrig %GPP toi is stim
                    cfg.latency = [0.050 0.150];
                    toi_epoch = ft_selectdata(cfg, data);
                end
                
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
                
                %invert GPPvsPO so "inhibition" is positive
                if ismember(plots(iii), 'GPPvsPO'); 
                    lcmv_contrast.avg.pow = (source_baseline.avg.pow - source_toi.avg.pow) ./ source_baseline.avg.pow;
                else
                    lcmv_contrast.avg.pow = (source_toi.avg.pow - source_baseline.avg.pow) ./ source_baseline.avg.pow;
                end
                
                %Interpolate contrast on MRI
                cfg = [];
                cfg.parameter = 'pow';
                cfg.interpmethod = 'nearest';
                
                source_int = ft_sourceinterpolate(cfg, lcmv_contrast, mri_resliced);
                
                %Find maxpower pos in POvsB
                if ismember(plots(iii), 'POvsB'); 
                    [M, I] = max(source_int.pow)
                    POmaxpos = source_int.pos(I,:);
                end

                figure('Position', [683,74,1095,1091])
                
                cfg = [];
                cfg.method = 'ortho';
                cfg.funparameter = 'pow';
                cfg.funcolormap = 'jet';
                
                %ft_sourceplot(cfg, source_int);
                
                %saveas(gcf, ['../Analysis Output/' name '_' sub_date.ID{i} '_' plots{iii} '.png'])
                close;

                %use POmaxpos for GPPvsPO contrast - else auto implied
                if ismember(plots(iii), 'GPPvsPO'); 

                    figure('Position', [683,74,1095,1091])

                    cfg = [];
                    cfg.method = 'ortho';
                    cfg.funparameter = 'pow';
                    cfg.funcolormap = 'jet';
                    
                    cfg.locationcoordinates = 'head';
                    cfg.location = POmaxpos;

                    ft_sourceplot(cfg, source_int);

                    saveas(gcf, ['../Analysis Output/' name '_' sub_date.ID{i} '_' plots{iii} '_POslice.png'])
                    close;
                end

            %For plots
            end

    %For bkgs
    end

%For subject
end

