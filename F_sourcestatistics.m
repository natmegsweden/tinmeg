
%outdir = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/'];

%% Load source reconstructions for conditions to cell array for all subjects

for i = 1:4%length(sub_date.ID);

inpath = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/'];

%Load conditions to compare
source_cond_PO60_75 = load([inpath 'PO60_70_stim_source.mat']); 
source_cond_PO60_75 = source_cond_PO60_75.stim_source;
PO60_75{i} = source_cond_PO60_75;

source_cond_PO60_90 = load([inpath 'PO60_90_stim_source.mat']); 
source_cond_PO60_90 = source_cond_PO60_90.stim_source;
PO60_90{i} = source_cond_PO60_90;

clear ('source_cond_PO60_75', 'source_cond_PO60_90');
end

%% Statistics over subjects

cfg=[];
cfg.dim         = PO60_75{1}.dim;
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.parameter   = 'pow';
cfg.correctm    = 'cluster';
cfg.numrandomization = 1000;
cfg.alpha       = 0.05; % note that this only implies single-sided testing
cfg.tail        = 0;

nsubj=numel(PO60_75);
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_sourcestatistics(cfg, PO60_75{:}, PO60_90{:});
