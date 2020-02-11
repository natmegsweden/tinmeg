addpath('C:\fieldtrip')
ft_defaults

%% 
trls = cleaned_downsampled_data.trialinfo==136 | cleaned_downsampled_data.trialinfo==132 ;

cfg =[];
cfg.trials = trls;
cfg.channel = 'megmag';
data_slct = ft_selectdata(cfg,  cleaned_downsampled_data);

tml = ft_timelockanalysis([],cleaned_downsampled_data);

tml_cmb = ft_combineplanar([], tml)

figure;
cfg = [];
cfg.layout = 'neuromag306planar';
ft_multiplotER(cfg, tml)

figure;
cfg = [];
cfg.layout = 'neuromag306cmb';
ft_multiplotER(cfg, tml_cmb)

figure;
cfg = [];
cfg.layout = 'neuromag306mag';
ft_multiplotER(cfg, tml_cmb)

%%
[val, idx] = max(max(tml_cmb.avg,[],2))
roi = tml_cmb.label{idx}

%% 

cfg =[];
cfg.trials = cleaned_downsampled_data.trialinfo==136;
data_136 = ft_selectdata(cfg,  cleaned_downsampled_data);
tml_136 = ft_timelockanalysis([], data_136);
tml_136_cmb = ft_combineplanar([], tml_136)

cfg.trials = cleaned_downsampled_data.trialinfo==132;
data_132 = ft_selectdata(cfg,  cleaned_downsampled_data);
tml_132 = ft_timelockanalysis([],data_132);
tml_132_cmb = ft_combineplanar([], tml_132)

figure;
cfg = [];
cfg.layout = 'neuromag306cmb';
ft_multiplotER(cfg, tml_cmb, tml_136_cmb, tml_132_cmb)


figure;
cfg = [];
cfg.layout = 'neuromag306mag';
ft_multiplotER(cfg, tml_cmb, tml_136_cmb, tml_132_cmb)

%% 
cfg = []
cfg.parameter = 'avg'
cfg.operation = 'x1-x2'
diff = ft_math(cfg, tml_136_cmb, tml_132_cmb)

figure;
cfg = [];
cfg.layout = 'neuromag306mag';
ft_multiplotER(cfg, diff,tml_136_cmb, tml_132_cmb)

%% 
cfg =[];
cfg.latency = [0.05 0.15];
cfg.channel = idx;

chan136 = ft_selectdata(cfg, tml_136_cmb)
chan132 = ft_selectdata(cfg, tml_132_cmb)

% Summaries
[valmax136, idxmax136] = max(chan136.avg);
chan136.time(idxmax136)
[valmin136, idxmin136] = min(chan136.avg);
chan136.time(idxmin136)

pk2pk136 = valmax136-valmin136;


[valmax132, idxmax132] = max(chan132.avg);
chan132.time(idxmax132)
[valmin136, idxmin136] = min(chan136.avg);
chan136.time(idxmin136)

pk2pk136 = valmax136-valmin136;



