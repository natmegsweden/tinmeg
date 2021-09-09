
%% Source statistics and analysis

brainnetome = ft_read_atlas('../../fieldtrip-20210311/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii', 'unit', 'mm');

%% Load source reconstructions for conditions to cell array for all subjects

for i = 1:4%length(sub_date.ID);

inpath = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/'];
    
    %PO60
    for ii = 1:length(cond.PO60label);

    %Load conditions to compare
    base = load([inpath cond.PO60label{ii} '_base_source.mat']);
    base = base.base_source;

    stim = load([inpath cond.PO60label{ii} '_stim_source.mat']);
    stim = stim.stim_source;
    
    %Overwrite position with template coordinates
    base.pos = template_grid.pos;
    stim.pos = template_grid.pos;
    
    %Write to structure
    PO60_base{i,ii} = base;
    PO60_stim{i,ii} = stim;

    end

end

clear ('base', 'stim');
%% Statistics over subjects

load standard_mri;
template_mri = mri;

for i = 1:length(PO60_stim(1,:));
    
    %PO60_stim{1,i};
    
    cfg=[];
    cfg.dim         = PO60_stim{1,i}.dim;
    cfg.method      = 'montecarlo';
    cfg.statistic   = 'ft_statfun_depsamplesT';
    cfg.parameter   = 'pow';
    cfg.correctm    = 'cluster';
    cfg.numrandomization = 1000;
    cfg.alpha       = 0.05; % note that this only implies single-sided testing
    cfg.tail        = 0;

    nsubj=numel(PO60_stim(:,i));
    cfg.design(1,:) = [1:nsubj 1:nsubj];
    cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
    cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
    cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

    stat = ft_sourcestatistics(cfg, PO60_stim{:,i}, PO60_base{:,i});
    
    cfg = [];
    cfg.parameter    = 'stat';
    cfg.interpmethod = 'nearest';
    cfg.downsample = 5;
    int_stat_ds = ft_sourceinterpolate(cfg, stat, template_mri);

    save(['../mat_data/' cond.PO60label{i} '_int_ds_stat.mat'], 'int_stat_ds', '-v7.3');
    
end

%%

cfg = [];
cfg.method          = 'slice';
cfg.funparameter    = 'stat';
cfg.funcolormap     = 'jet';

cfg.maskparameter = cfg.funparameter;

cfg.funcolorlim = [-5 0];

cfg.opacitylim = [-5 0];
cfg.opacitymap = 'rampup';

ft_sourceplot(cfg, stat_temp_int);


