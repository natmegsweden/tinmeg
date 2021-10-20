%% Source statistics and analysis

%load atlas
brainnetome = ft_read_atlas('../../fieldtrip-20210311/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii', 'unit', 'mm');
%aal_atlas = ft_read_atlas('../../fieldtrip-20210311/template/atlas/aal/ROI_MNI_V4.nii', 'unit', 'mm');

%Interpolate atlas to template_grid
cfg = [];
cfg.parameter = 'tissue';
cfg.interpmethod = 'nearest';
atlas_grid = ft_sourceinterpolate(cfg, brainnetome, template_grid);

atlas_grid = ft_checkdata(atlas_grid, 'datatype', 'source');
atlas_grid.inside = template_grid.inside;

%Brainnetome-atlas labels (Brodmann 41&42 + TETE1.0&TE1.2)
labels = [71 72 73 74];

%Plot Atlas-ROI
cfg.parameter = 'tissue';
cfg.interpmethod = 'nearest';
atlas_grid_int = ft_sourceinterpolate(cfg, atlas_grid, template_mri);

dummymask = (ismember(atlas_grid_int.tissue, labels));
atlas_grid_int.dummymask = dummymask;

cfg = [];
cfg.method          = 'ortho';
cfg.funparameter    = 'tissue';
cfg.anaparameter    = 'anatomy';
cfg.maskparameter   = 'dummymask';
cfg.nslices = 25;
cfg.position        = [700 300 950 950];
ft_sourceplot(cfg, atlas_grid_int);

%% Plot difference map in TOI - WIP

    powdiff_sub = load([inpath cond.GP60label{3} '_powdiff.mat']);
    powdiff_sub = powdiff_sub.powdiff_sub;
    powdiff_sub.pos = template_grid.pos;
    
    %POWERDIFF - Works but not informative..
    cfg.parameter = 'pow';
    cfg.interpmethod = 'nearest';
    pow_diff_int = ft_sourceinterpolate(cfg, powdiff_sub, template_mri);
    
    pow_diff_int.dummymask = dummymask;
    
    cfg = [];
    cfg.method          = 'slice';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';
    cfg.maskparameter   = 'dummymask';
    cfg.position        = [700 300 950 950];
    ft_sourceplot(cfg, pow_diff_int);

%% Plot virtual channel for ROI

%virtchans = struct;
%load virtchans

%Condition
for iii = 1%length(conditions);

    %Trial
    for i = 5%1:length(cond.([conditions{iii} 'label']));

        for ii = 1:length(sub_date.ID)

        inpath = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{ii} '/'];

        stim = load([inpath cond.PO60label{i} '_stim_source.mat']);
        stim = stim.stim_source; 
        base = load([inpath cond.PO60label{i} '_base_source.mat']);
        base = base.base_source;

        %Parcellation
        stim.tissue = atlas_grid.tissue;
        stim.tissuelabel = atlas_grid.tissuelabel;
        stim.pos = template_grid.pos;

        base.tissue = atlas_grid.tissue;
        base.tissuelabel = atlas_grid.tissuelabel;
        base.pos = template_grid.pos;

        cfg = [];
        cfg.parcellation = 'tissue';
        cfg.method = 'eig';

        sub_par_stim = ft_sourceparcellate(cfg, stim, atlas_grid);
        sub_par_base = ft_sourceparcellate(cfg, base, atlas_grid);

        %Virtual-channel for area 71:74 in brainnetome
        virtchans.(cond.([conditions{iii} 'label']){i})(ii,1:61) = mean(sub_par_stim.mom(71:74, :)) - (mean(mean(sub_par_base.mom(71:74, :))));

        clear stim base;
        
        %Subject
        end

    %Trials
    end
    
%Conditions
end

%Condition
for iii = 3%length(conditions);

    %Trial
    for i = 3%1:length(cond.([conditions{iii} 'label']));

        for ii = 1:length(sub_date.ID)

        inpath = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{ii} '/'];

        stim = load([inpath cond.GP60label{i} '_stim_source.mat']);
        stim = stim.stim_source; 
        base = load([inpath cond.GP60label{i} '_base_source.mat']);
        base = base.base_source;

        %Parcellation
        stim.tissue = atlas_grid.tissue;
        stim.tissuelabel = atlas_grid.tissuelabel;
        stim.pos = template_grid.pos;

        base.tissue = atlas_grid.tissue;
        base.tissuelabel = atlas_grid.tissuelabel;
        base.pos = template_grid.pos;

        cfg = [];
        cfg.parcellation = 'tissue';
        cfg.method = 'eig';

        sub_par_stim = ft_sourceparcellate(cfg, stim, atlas_grid);
        sub_par_base = ft_sourceparcellate(cfg, base, atlas_grid);

        %Virtual-channel for area 71:74 in brainnetome
        virtchans.(cond.([conditions{iii} 'label']){i})(ii,1:61) = mean(sub_par_stim.mom(71:74, :)) - (mean(mean(sub_par_base.mom(71:74, :))));

        clear stim base;
        
        %Subject
        end

    %Trials
    end
    
%Conditions
end

%save(['../mat_data/source_reconstruction/'  'virtchans_i120.mat'], 'virtchans', '-v7.3');

figure; hold on
for p = 1:22
   
    subplot(2,1,1);
    xlim([1 61]);
    ylim([-4 3]);
    title('Pulse only control');
    plot(virtchans.PO60_90(p,:), 'Color', [1 0 0 0.25]); hold on
    plot(mean(virtchans.PO60_90(:,:)), 'Color', [1 0 0 1], 'LineWidth', 1);
    
    x = gca;
    x.XTick = [1:10:61];
    x.XTickLabel = [0:50:300];
    x.XMinorTick = 'on';
    x.FontSize = 12;
    
    subplot(2,1,2);
    xlim([1 61]);
    ylim([-4 3]);
    title('ISI 120');
    plot(virtchans.GP60_i120(p,:), 'Color', [0 0 1 0.25]); hold on;
    plot(mean(virtchans.GP60_i120(:,:)), 'Color', [0 0 1 1], 'LineWidth', 1);
    
    x = gca;
    x.XTick = [1:10:61];
    x.XTickLabel = [0:50:300];
    x.XMinorTick = 'on';
    x.FontSize = 12;
    
end
    
%%

all_stim{ii} = stim;
all_base{ii} = base;

    %Calculate and save power difference for grand average
    cfg = [];
    cfg.keepindividual = 'no';

    stim_gravg = ft_sourcegrandaverage(cfg, all_stim{:});
    base_gravg = ft_sourcegrandaverage(cfg, all_base{:});
    
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = 'x1 - x2';
    pow_diff = ft_math(cfg, stim_gravg, base_gravg);
    
    clear('stim_gravg', 'base_gravg');
    
    %save(['../mat_data/stats/' cond.PO60label{i} '_pow_diff.mat'], 'pow_diff', '-v7.3');
    
    
    %Sourcestatistic between stim and baseline window
    cfg=[];
    cfg.dim         = all_stim{1}.dim;
    cfg.method      = 'montecarlo';
    cfg.statistic   = 'ft_statfun_depsamplesT';
    cfg.parameter   = 'pow';
    cfg.correctm    = 'cluster';
    cfg.numrandomization = 'all';
    cfg.alpha       = 0.05; %corrected below
    cfg.tail        = 0; %two-sided
    cfg.correcttail = 'alpha'; %http://bit.ly/2YQ1Hm8 

    nsubj=numel(all_stim);
    cfg.design(1,:) = [1:nsubj 1:nsubj];
    cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
    cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
    cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

    stat = ft_sourcestatistics(cfg, all_stim{:}, all_base{:});
    
    %save(['../mat_data/stats/' cond.PO60label{i} '_stat_stimVSbase.mat'], 'stat', '-v7.3');
    
    clear ('all_base', 'all_stim');
    
    %Interpolate power difference on template MRI
    cfg = [];
    cfg.parameter    = 'pow';
    cfg.interpmethod = 'nearest';
    int_powdiff = ft_sourceinterpolate(cfg, pow_diff, template_mri);

    %save(['../mat_data/stats/' cond.PO60label{i} '_int_powdiff.mat'], 'int_powdiff', '-v7.3');
    
    %clear ('int_powdiff');
    

 
%ORTHO NOMASK
cfg = [];
cfg.method          = 'ortho';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';

%cfg.maskparameter = 'statmask';

% cfg.funcolorlim = [0 20];
% cfg.opacitylim = [0 20];
cfg.opacitymap = 'rampup';

cfg.position = [700 300 950 950];
ft_sourceplot(cfg, int_powdiff);
