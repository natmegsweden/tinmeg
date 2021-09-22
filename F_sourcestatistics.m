
%% Source statistics and analysis

brainnetome = ft_read_atlas('../../fieldtrip-20210311/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii', 'unit', 'mm');

load standard_mri;
template_mri = mri;

%% Load subjects source reconstructions for conditions and run ft_sourcestatistic
    
    %PO60
    %cond.PO60label is the names of 6 conditions i.e. 'PO60_70', 'PO60_75' and so on..
    
for i = 1%:length(cond.PO60label);

    for ii = 1%:4%length(sub_date.ID) - load one subject only for test
    
    inpath = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{ii} '/']; %sub_date is table of subject IDs
        
    %Load conditions to compare
    base = load([inpath cond.PO60label{i} '_base_source.mat']); %baseline from ft_sourceanalysis
    base = base.base_source;

    stim = load([inpath cond.PO60label{i} '_stim_source.mat']); %stim response window from ft_sourceanalysis
    stim = stim.stim_source;
    
    %Overwrite position with template coordinates
    base.pos = template_grid.pos;
    stim.pos = template_grid.pos;
    
    %Collect all subjects in structure
    all_stim{ii} = stim
    all_base{ii} = base
    
    clear('stim', 'base');
    
    end
    
    %Calculate and save power difference for grand average
    cfg = [];
    cfg.keepindividual = 'no';

    stim_gravg = ft_sourcegrandaverage(cfg, all_stim{:});
    base_gravg = ft_sourcegrandaverage(cfg, all_base{:});
    
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = 'x1 - x2'
    pow_diff = ft_math(cfg, stim_gravg, base_gravg);
    
    clear('stim_gravg', 'base_gravg');
    
    %save(['../mat_data/stats/' cond.PO60label{i} '_pow_diff.mat'], 'pow_diff', '-v7.3');
    
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
    
    cfg = [];
    cfg.parameter    = 'pow';
    cfg.interpmethod = 'nearest';
    int_powdiff = ft_sourceinterpolate(cfg, pow_diff, template_mri);

    %save(['../mat_data/stats/' cond.PO60label{i} '_int_powdiff.mat'], 'int_powdiff', '-v7.3');
    
    %clear ('int_powdiff');
    
end
 
%% Plot source power difference, n=4, PO60

for i = 1:numel(cond.PO60label);

load(['../mat_data/stats/'  cond.PO60label{i} '_stat_stimVSbase.mat']);

%interpolate %.mask
cfg = [];
cfg.parameter    = 'mask';
cfg.interpmethod = 'nearest';
int_statmask = ft_sourceinterpolate(cfg, stat, template_mri);

clear stat;

save(['../mat_data/stats/' cond.PO60label{i} '_int_statmask.mat'], 'int_statmask', '-v7.3');

load(['../mat_data/stats/'  cond.PO60label{i} '_int_powdiff.mat']);

%Reshape .mask into cell array in interpolated power difference
statmask = reshape(int_statmask.mask, int_powdiff.dim);
int_powdiff.statmask = statmask;

clear int_statmask;

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

saveas(gcf, ['../Analysis Output/source_recon_test/' cond.PO60label{i} '_ortho_nomask.svg']);

close


%ORTHO MASK
cfg = [];
cfg.method          = 'ortho';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';

cfg.maskparameter = 'statmask';

% cfg.funcolorlim = [0 20];
% cfg.opacitylim = [0 20];
cfg.opacitymap = 'rampup';

cfg.position = [700 300 950 950];
ft_sourceplot(cfg, int_powdiff);

saveas(gcf, ['../Analysis Output/source_recon_test/' cond.PO60label{i} '_ortho_mask.svg']);

close


%SLICE NOMASK
cfg = [];
cfg.method          = 'slice';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';

%cfg.maskparameter = 'statmask';

% cfg.funcolorlim = [0 20];
% cfg.opacitylim = [0 20];
cfg.opacitymap = 'rampup';

cfg.position = [700 300 950 950];
ft_sourceplot(cfg, int_powdiff);

saveas(gcf, ['../Analysis Output/source_recon_test/' cond.PO60label{i} '_slice_nomask.svg']);

close


%SLICE MASK
cfg = [];
cfg.method          = 'slice';
cfg.funparameter    = 'pow';
cfg.funcolormap     = 'jet';

cfg.maskparameter = 'statmask';

% cfg.funcolorlim = [0 20];
% cfg.opacitylim = [0 20];
cfg.opacitymap = 'rampup';

cfg.position = [700 300 950 950];
ft_sourceplot(cfg, int_powdiff);

saveas(gcf, ['../Analysis Output/source_recon_test/' cond.PO60label{i} '_slice_mask.svg']);

close

end
