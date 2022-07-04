
load('../mat_data/timelockeds/all_cmb_avg.mat');

% need to run ft_timelockgrandaverage to keep label field
% see https://www.fieldtriptoolbox.org/tutorial/eventrelatedstatistics/#permutation-test-based-on-cluster-statistics

cfg = [];
cmb_timelockeds = ft_combineplanar(cfg, timelockeds);

cfg = [];
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';
cfg.layout = 'neuromag306cmb_neighb.mat';
cfg.feedback    = 'yes';                              % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, cmb_timelockeds); % define neighbouring channels

ft_timelockstatistics()