

% For every subject
for i = 1:length(sub_date.ID(1:4)) %test with first 4
    
%load and append MEG data
%infiles = find_files(['../mat_data/preprocessing/' 'ID' sub_date.ID{i}], 'ds_clean');
%fname = eval(sub_date.ID{i}) ??

%Load in the stupidest way possible
load(['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' 'GO_ds_clean.mat']);
GOmat = cleaned4mat;

load(['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' 'PO60_ds_clean.mat']);
PO60mat = cleaned4mat;

load(['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' 'PO70_ds_clean.mat']);
PO70mat = cleaned4mat;

load(['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' 'GP60_ds_clean.mat']);
GP60mat = cleaned4mat;

load(['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' 'GP70_ds_clean.mat']);
GP70mat = cleaned4mat;

cfg = [];
appended = ft_appenddata(cfg, PO60mat, PO70mat, GP60mat, GP70mat, GOmat);

end



cfg              = [];
cfg.method       = 'fastica';
cfg.numcomponent = 40;
cfg.channel      = 'MEG';
comp             = ft_componentanalysis(cfg, cleaned4mat);

%% paj

cfg = [];
cfg.layout      = 'neuromag306all.lay';
cfg.marker      = 'off';
cfg.component   = [1:20];               
figure; ft_topoplotIC(cfg, comp);

% Time-series view (split in two figures)
cfg = [];
cfg.viewmode    = 'component';
cfg.layout      = 'neuromag306all.lay';
cfg.blocksize   = 10;
cfg.channel     = [1:20];
ft_databrowser(cfg, comp)
%% Semiauto ECG ICA

%Select QRS threshold and time window
cfg = [];
cfg.continuous            = 'no';
cfg.artfctdef.ecg.pretim  = 0.075; %0.3
cfg.artfctdef.ecg.psttim  = 0.1; %0.3
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
cfg.channel    = 'MEG*'; %??
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

%% Find ECG components
% Define cutoff
cutoff = 0.5;           % Between 0-1 (analogue to a correlation coefficient)

% compute a frequency decomposition of all components and the ECG
cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.foilim     = [0 30];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq = ft_freqanalysis(cfg, comp_ecg);

% compute coherence between all components and the ECG
cfg = [];
cfg.channelcmb = {'all' 'ECG'};
cfg.method     = 'coh';
fdcomp = ft_connectivityanalysis(cfg, freq);

% Find ECG components
maxcoh = max(fdcomp.cohspctrm, [], 2);
ecg_comp_idx = find(maxcoh > cutoff);

%%

figure;
subplot(3,1,1); 
plot(fdcomp.freq, abs(fdcomp.cohspctrm)); hold on
plot([min(fdcomp.freq),max(fdcomp.freq)],[cutoff, cutoff], 'k--')
title('ECG'); xlabel('freq'); ylabel('coh');

subplot(3,1,2); imagesc(abs(fdcomp.cohspctrm));
xlabel('freq'); ylabel('comp');
subplot(3,1,3);
maxcoh = max(fdcomp.cohspctrm, [], 2);
foo = find(~(maxcoh > cutoff));
bp = bar(1:length(maxcoh), diag(maxcoh), 'stacked');
set(bp(foo),'facecolor','w'); set(bp(ecg_comp_idx),'facecolor','r')
axis([0.5, length(maxcoh)+0.5, 0, 1]); xlabel('comp'); ylabel('coh');

% View marked component(s)
cfg = [];
cfg.channel     = ecg_comp_idx; % components to be plotted
cfg.viewmode    = 'component';
cfg.layout      = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
ft_databrowser(cfg, comp)

%% reject component

%ecg_comp_idx:
reject_comp = [1, 3, 7, 11];    % Write the index of the components you want to remove

% Remove components
cfg = [];
cfg.component   = reject_comp; %ecg_comp_idx
cfg.channel     = 'MEG';
cfg.updatesens  = 'no'; %
icacleaned_downsampled_data = ft_rejectcomponent(cfg, comp, cleaned_downsampled_data);