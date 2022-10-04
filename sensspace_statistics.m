%% Load subjects timlocked data

%Cells for conditions of interest

PO60_90 = cell(1,1);
PO70_90 = cell(1,1);

GP60_i240 = cell(1,1);
GP70_i240 = cell(1,1);

% NB - MOVE ft_combineplanar to sensspace_processing an load directly in
% below loops. Save cell arrays somewhere.

% %Load each subject to the four conditions and collect to cell array
% %PO60_90
% for i = 1:numel(sub_date.ID)
%     temp = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/PO60_90_tlks.mat']);
%     temp = temp.timelockeds;
% 
%     %%%%%%%%%%%%%%%%%% temporary workaround
%     cfg = [];
%     temp = ft_combineplanar(cfg, temp);
%     %%%%%%%%%%%%%%%%%% temporary workaround
% 
%     PO60_90{1,i} = temp;
%     clear temp
% end
% 
% %GP60_i240
% for i = 1:numel(sub_date.ID)
%     temp = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/GP60_i240_tlks.mat']);
%     temp = temp.timelockeds;
% 
%     %%%%%%%%%%%%%%%%%% temporary workaround
%     cfg = [];
%     temp = ft_combineplanar(cfg, temp);
%     %%%%%%%%%%%%%%%%%% temporary workaround
% 
%     GP60_i240{1,i} = temp;
%     clear temp
% end

%PO70_90
for i = 1:numel(sub_date.ID)
    temp = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/PO70_90_tlks.mat']);
    temp = temp.timelockeds;

    %%%%%%%%%%%%%%%%%% temporary workaround
    cfg = [];
    temp = ft_combineplanar(cfg, temp);
    %%%%%%%%%%%%%%%%%% temporary workaround

    PO70_90{1,i} = temp;
    clear temp
end

%GP70_i240
for i = 1:numel(sub_date.ID)
    temp = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/GP70_i240_tlks.mat']);
    temp = temp.timelockeds;

    %%%%%%%%%%%%%%%%%% temporary workaround
    cfg = [];
    temp = ft_combineplanar(cfg, temp);
    %%%%%%%%%%%%%%%%%% temporary workaround

    GP70_i240{1,i} = temp;
    clear temp
end

%Prepare struct describing neighbouring channels
cfg = [];
cfg.method = 'distance';
cfg.channel = 'MEG';

neighbours = ft_prepare_neighbours(cfg, PO60_90{1,1});


%% Cluster based permutation statistic
cfg         = [];
cfg.channel = {'MEG'};
cfg.latency = [-0.5 0.32];

cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
%cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

Nsubj  = numel(sub_date.ID);
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

[stat_allchan] = ft_timelockstatistics(cfg, PO60_90{:}, GP60_i240{:});

cfg.channel = top_chan60;
cfg.latency = [0.05 0.15];

[stat_soichan] = ft_timelockstatistics(cfg, PO60_90{:}, GP60_i240{:});

stat_allchan.description = "Testing PO60_90 vs GP60_i240 for all MEG sensors over full trial window (-500 to 320ms)";
stat_soichan.description = "Testing PO60_90 vs GP60_i240 for grads of interest from sensspace_analysis in TOI for N1 only";

% >2GB file size :O
%save('../mat_data/stats/clustertest_allchan.mat', 'stat_allchan', '-v7.3')
%save('../mat_data/stats/clustertest_soichan.mat', 'stat_soichan', '-v7.3')

p_allchan = stat_allchan.prob;
p_soichan = stat_soichan.prob;

save('../mat_data/stats/clustertest_p_allchan.mat', 'p_allchan');
save('../mat_data/stats/clustertest_p_soichan.mat', 'p_soichan');

%% 70 carrier Cluster based permutation statistic
cfg         = [];
cfg.channel = {'MEG'};
cfg.latency = [-0.5 0.32];

cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
%cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

Nsubj  = numel(sub_date.ID);
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

[stat_allchan] = ft_timelockstatistics(cfg, PO70_90{:}, GP70_i240{:});

cfg.channel = top_chan70;
cfg.latency = [0.05 0.15];

[stat_soichan] = ft_timelockstatistics(cfg, PO70_90{:}, GP70_i240{:});

stat_allchan.description = "Testing PO70_90 vs GP70_i240 for all MEG sensors over full trial window (-500 to 320ms)";
stat_soichan.description = "Testing PO70_90 vs GP70_i240 for grads of interest from sensspace_analysis in TOI for N1 only";

% >2GB file size :O
%save('../mat_data/stats/clustertest_allchan.mat', 'stat_allchan', '-v7.3')
%save('../mat_data/stats/clustertest_soichan.mat', 'stat_soichan', '-v7.3')

p_allchan70 = stat_allchan.prob;
p_soichan70 = stat_soichan.prob;

save('../mat_data/stats/clustertest_p_allchan70.mat', 'p_allchan70');
save('../mat_data/stats/clustertest_p_soichan70.mat', 'p_soichan70');
