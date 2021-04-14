%% Load downsampled data and clean with ft_rejectvisual

for i = 1:height(sub_date)

fname = ['ID' char(fpaths(i,1)) '_60PO_ds_clean' '.mat'];
fpath = ['../mat_data/' fname];
    
%check if file exist
if exist(fpath, 'file')
    warning(['Output exist for subject ' char(fpaths(i,1))])
    continue
end

fname = ['ID' char(fpaths(i,1)) '_60PO_ds' '.mat']
fpath = ['../mat_data/' fname]

cleaned4mat = load(fpath);
cleaned4mat = cleaned4mat.res4mat_ds

cfg = [];
cfg.method = 'summary';
cfg.keepchannel = 'yes';
cfg.channel = 'MEGMAG';
cfg.layout = 'neuromag306eeg_1005_natmeg.lay';

cleaned4mat = ft_rejectvisual(cfg,cleaned4mat);

cfg.channel = 'MEGGRAD';

cleaned4mat = ft_rejectvisual(cfg, cleaned4mat);

fname = ['ID' char(fpaths(i,1)) '_60PO_ds_clean' '.mat'];
fpath = ['../mat_data/' fname];

save(fpath, 'cleaned4mat');
end

%ft_badchannel + ft_badsegment, or log from ft_rejectvisual?