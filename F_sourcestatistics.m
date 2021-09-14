
%% Source statistics and analysis

brainnetome = ft_read_atlas('../../fieldtrip-20210311/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii', 'unit', 'mm');

load standard_mri;
template_mri = mri;

%% Load subjects source reconstructions for conditions and run ft_sourcestatistic
    
    %PO60
for i = 1:length(cond.PO60label);

    for ii = 1:4%length(sub_date.ID)
        
    inpath = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{ii} '/'];
        
    %Load conditions to compare
    base = load([inpath cond.PO60label{i} '_base_source.mat']);
    base = base.base_source;

    stim = load([inpath cond.PO60label{i} '_stim_source.mat']);
    stim = stim.stim_source;
    
    %Overwrite position with template coordinates
    base.pos = template_grid.pos;
    stim.pos = template_grid.pos;
    
    all_stim{ii} = stim
    all_base{ii} = base
    
    end

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
    
    save(['../mat_data/stats/' cond.PO60label{i} '_stat_stimVSbase.mat'], 'stat', '-v7.3');
    
    cfg = [];
    cfg.parameter    = 'stat';
    cfg.interpmethod = 'nearest';
    int_stat = ft_sourceinterpolate(cfg, stat, template_mri);

    save(['../mat_data/stats/' cond.PO60label{i} '_int_stat.mat'], 'int_stat', '-v7.3');
    
end



clear ('base', 'stim');

%%

cfg = [];
cfg.method          = 'slice';
cfg.funparameter    = 'stat';
cfg.funcolormap     = 'jet';

cfg.maskparameter = cfg.funparameter;

cfg.funcolorlim = [0 250];

cfg.opacitylim = [0 250];
cfg.opacitymap = 'rampup';

ft_sourceplot(cfg, int_stat);


