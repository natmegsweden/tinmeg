
%%

load('../mat_data/preprocessing/ID0697/PO60_ds_clean.mat');

%%
cfg            = [];
cfg.method     = 'runica';
comp           = ft_componentanalysis(cfg, cleaned4mat);

%% paj

cfg = [];
cfg.layout      = 'neuromag306all.lay';
cfg.marker      = 'off';
cfg.component   = [1:20];               
figure; ft_topoplotIC(cfg, comp);

%% Semiauto ECG ICA

%Select QRS threshold and time window
cfg = [];
cfg.continuous            = 'no';
cfg.artfctdef.ecg.pretim  = 0.25;
cfg.artfctdef.ecg.psttim  = 0.50;
cfg.channel               = 'ECG';
cfg.artfctdef.ecg.inspect = 'ECG';
[~, artifact] = ft_artifact_ecg(cfg, cleaned4mat);

% Make artifact epochs
cfg = [];
cfg.dftfilter  = 'yes';
cfg.demean     = 'yes';
cfg.trl        = [artifact zeros(size(artifact,1), 1)];
temp = ft_redefinetrial(cfg, cleaned4mat);

% Re-arrange data
cfg.channel    = 'MEG*';
data_ecg = ft_selectdata(cfg, temp);

cfg.channel    = 'ECG';
ecg = ft_selectdata(cfg, temp);

ecg.channel{:} = 'ECG';         % renaming for bookkeeping

% Remove residual linenoise in electric channel.
cfg = [];
cfg.dftfilter       = 'yes';
cfg.dftfreq         = [50, 100, 150];
ecg = ft_preprocessing(cfg, ecg);

% decompose the ECG-locked datasegments (using the unmixing matrix from comp)
cfg = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_ecg = ft_componentanalysis(cfg, data_ecg);

% append the ecg channel to the data structure;
comp_ecg = ft_appenddata([], ecg, comp_ecg);

% average the components timelocked to the QRS-complex
cfg = [];
timelock = ft_timelockanalysis(cfg, comp_ecg);

% Plot
figure
subplot(2,1,1); plot(timelock.time, timelock.avg(1,:)); title('ECG')
subplot(2,1,2); plot(timelock.time, timelock.avg(2:end,:));  title('ICA comp')

