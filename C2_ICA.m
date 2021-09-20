

% For every subject

inpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/']

outdir = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/']

%Check if subject dir exist, create/define
if ~exist(outdir, 'file');
    mkdir(outdir);
end

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

%Append data
cfg = [];
cfg.keepsampleinfo = 'no';
appended = ft_appenddata(cfg, PO60mat, PO70mat, GP60mat, GP70mat, GOmat);

%Identify components and save
cfg              = [];
cfg.method       = 'fastica';
cfg.numcomponent = 40;
cfg.channel      = 'MEG';
comp             = ft_componentanalysis(cfg, appended);

save([outdir 'components.mat'], 'comp', '-v7.3');

% %Plot
% cfg = [];
% cfg.layout      = 'neuromag306all.lay';
% cfg.marker      = 'off';
% cfg.component   = [1:20];               
% figure; ft_topoplotIC(cfg, comp);
% 
% % Time-series view (split in two figures)
% cfg = [];
% cfg.viewmode    = 'component';
% cfg.layout      = 'neuromag306all.lay';
% cfg.blocksize   = 10;
% cfg.channel     = [1:20];
% ft_databrowser(cfg, comp)

%Select QRS threshold and time window, requires inspection
cfg = [];
cfg.continuous            = 'no';
cfg.artfctdef.ecg.pretim  = 0.3; %0.3
cfg.artfctdef.ecg.psttim  = 0.3; %0.3
cfg.artfctdef.ecg.cutoff  = 3; %Default 3
cfg.channel               = 'ECG';
cfg.artfctdef.ecg.inspect = 'ECG';
[~, artifact] = ft_artifact_ecg(cfg, appended);

% Make artifact epochs
cfg = [];
cfg.dftfilter  = 'yes';
cfg.demean     = 'yes';
cfg.trl        = [artifact zeros(size(artifact,1), 1)];
temp = ft_redefinetrial(cfg, appended);

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

save([outdir 'components_ecg.mat'], 'comp_ecg', '-v7.3');

% average the components timelocked to the QRS-complex to inspect
% cfg = [];
% timelock = ft_timelockanalysis(cfg, comp_ecg);

% % Plot
% figure
% subplot(2,1,1); plot(timelock.time, timelock.avg(1,:)); title('ECG')
% subplot(2,1,2); plot(timelock.time, timelock.avg(2:end,:));  title('ICA comp')



% Find ECG components

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

%Create and save plots of components
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

saveas(gcf, [outdir 'comp_figure.pdf'])

close(1)

% % View marked component(s)
% cfg = [];
% cfg.channel     = ecg_comp_idx; % components to be plotted
% cfg.viewmode    = 'component';
% cfg.layout      = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
% ft_databrowser(cfg, comp)


% reject component

% Remove components
cfg = [];
cfg.component   = ecg_comp_idx; %index of components to remove
cfg.channel     = 'MEG';
cfg.updatesens  = 'no'; %
appended_ICA = ft_rejectcomponent(cfg, comp, appended);


%Split conditions and save

%NB, order from ft_appenddata:
%PO60mat, PO70mat, GP60mat, GP70mat, GOmat

%PO60
first = 1;
last = length(PO60mat.trial);

cfg = [];
cfg.trials = (first:last);
PO60ica = ft_selectdata(cfg, appended_ICA);
PO60ica.sampleinfo = PO60mat.sampleinfo;

save([outdir 'PO60ica.mat'], 'PO60ica');

%PO70
first = last+1;
last = first+length(PO70mat.trial)-1;

cfg = [];
cfg.trials = (first:last);
PO70ica = ft_selectdata(cfg, appended_ICA);
PO70ica.sampleinfo = PO70mat.sampleinfo;

save([outdir 'PO70ica.mat'], 'PO70ica');

%GP60
first = last+1;
last = first+length(GP60mat.trial)-1;

cfg = [];
cfg.trials = (first:last);
GP60ica = ft_selectdata(cfg, appended_ICA);
GP60ica.sampleinfo = GP60mat.sampleinfo;

save([outdir 'GP60ica.mat'], 'GP60ica');

%GP70
first = last+1;
last = first+length(GP70mat.trial)-1;

cfg = [];
cfg.trials = (first:last);
GP70ica = ft_selectdata(cfg, appended_ICA);
GP70ica.sampleinfo = GP70mat.sampleinfo;

save([outdir 'GP70ica.mat'], 'GP70ica');

%GO
first = last+1;
last = first+length(GOmat.trial)-1;

cfg = [];
cfg.trials = (first:last);
GOica = ft_selectdata(cfg, appended_ICA);
GOica.sampleinfo = GOmat.sampleinfo;

save([outdir 'GOica.mat'], 'GOica');

clear('PO60ica', 'PO70ica', 'GP60ica', 'GP70ica', 'GOica', 'appended_ICA', 'appended', 'temp', 'data_ecg', 'comp_ecg', 'freq', 'fd_comp', 'fdcomp', 'comp', 'PO60mat', 'PO70mat', 'GP60mat', 'GP70mat', 'GOmat', 'cleaned4mat', 'ecg');
