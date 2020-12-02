addpath('/home/mikkel/fieldtrip/fieldtrip/')
ft_defaults

%% Load data



%% Stats
% Make the design matric
ivar = [60 65 70 75]; % The predictor
uvar = [1 1 1 1];     % Unit of ubervation (i.e. within subject identifier)

design = [ivar', uvar'];

% Use stats stat
cfg = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.statistic           = 'ft_statfun_depsamplesregrT';
cfg.alpha               = 0.025;
cfg.tail                = 0;
cfg.clusterstatistic    = 'maxsum';
cfg.clusteralpha        = 0.05;
cfg.clustertail         = 0;            % = two-tailed
cfg.numrandomization    = 1000;
cfg.ivar                = 1;
cfg.uvar                = 2;            % Specific for dependent tests
cfg.design              = design;

stats = ft_timelockstatistics(cfg, data)



