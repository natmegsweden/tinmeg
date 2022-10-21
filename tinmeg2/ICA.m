% For every subject

inpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/']
outdir = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/']

%Check if subject dir exist, create/define
if ~exist(outdir, 'file');
    mkdir(outdir);
end


%Load preprocessed data
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin0_bkg0_ds_clean.mat']);
tin0_bkg0 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin0_bkg3_ds_clean.mat']);
tin0_bkg3 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin0_bkg8_ds_clean.mat']);
tin0_bkg8 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin3_bkg0_ds_clean.mat']);
tin3_bkg0 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin3_bkg3_ds_clean.mat']);
tin3_bkg3 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin3_bkg8_ds_clean.mat']);
tin3_bkg8 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin8_bkg0_ds_clean.mat']);
tin8_bkg0 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin8_bkg3_ds_clean.mat']);
tin8_bkg3 = pre_tinmeg2_c;
load(['../mat_data/preprocessing/' 'ID' num2str(sub_date.ID{i}) '/' 'tin8_bkg8_ds_clean.mat']);
tin8_bkg8 = pre_tinmeg2_c;

clear tinmeg2_pre_c

%Append data
cfg = [];
cfg.keepsampleinfo = 'no';
appended = ft_appenddata(cfg, tin0_bkg0, tin0_bkg3, tin0_bkg8, tin3_bkg0, tin3_bkg3, tin3_bkg8, tin8_bkg0, tin8_bkg3, tin8_bkg8);

%Identify components and save
cfg              = [];
cfg.method       = 'fastica';
cfg.numcomponent = 40;
cfg.channel      = 'MEG';
comp             = ft_componentanalysis(cfg, appended);

save([outdir 'components_tinmeg2.mat'], 'comp', '-v7.3');

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

save([outdir 'components_ecg_tinmeg2.mat'], 'comp_ecg', '-v7.3');

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

saveas(gcf, [outdir 'comp_figure_tinmeg2.pdf'])

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
%tin0_bkg0, tin0_bkg3, tin0_bkg8, tin3_bkg0, tin3_bkg3, tin3_bkg8, tin8_bkg0, tin8_bkg3, tin8_bkg8

%tin0_bkg0
first = 1;
last = numel(tin0_bkg0.trial);

cfg = [];
cfg.trials = (first:last);
tin0_bkg0_ica = ft_selectdata(cfg, appended_ICA);
tin0_bkg0_ica.sampleinfo = tin0_bkg0.sampleinfo;
save([outdir 'tin0_bkg0_ica.mat'], 'tin0_bkg0_ica');

%tin0_bkg3
first = last+1;
last = first+numel(tin0_bkg3.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin0_bkg3_ica = ft_selectdata(cfg, appended_ICA);
tin0_bkg3_ica.sampleinfo = tin0_bkg3.sampleinfo;
save([outdir 'tin0_bkg3_ica.mat'], 'tin0_bkg3_ica');

%tin0_bkg8
first = last+1;
last = first+length(tin0_bkg8.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin0_bkg8_ica = ft_selectdata(cfg, appended_ICA);
tin0_bkg8_ica.sampleinfo = tin0_bkg8.sampleinfo;
save([outdir 'tin0_bkg8_ica.mat'], 'tin0_bkg8_ica');

%tin3_bkg0
first = last+1;
last = first+length(tin3_bkg0.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin3_bkg0_ica = ft_selectdata(cfg, appended_ICA);
tin3_bkg0_ica.sampleinfo = tin3_bkg0.sampleinfo;
save([outdir 'tin3_bkg0_ica.mat'], 'tin3_bkg0_ica');

%tin3_bkg3
first = last+1;
last = first+length(tin3_bkg3.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin3_bkg3_ica = ft_selectdata(cfg, appended_ICA);
tin3_bkg3_ica.sampleinfo = tin3_bkg3.sampleinfo;
save([outdir 'tin3_bkg3_ica.mat'], 'tin3_bkg3_ica');

%tin3_bkg8
first = last+1;
last = first+length(tin3_bkg8.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin3_bkg8_ica = ft_selectdata(cfg, appended_ICA);
tin3_bkg8_ica.sampleinfo = tin3_bkg8.sampleinfo;
save([outdir 'tin3_bkg8_ica.mat'], 'tin3_bkg8_ica');

%tin8_bkg0
first = last+1;
last = first+length(tin8_bkg0.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin8_bkg0_ica = ft_selectdata(cfg, appended_ICA);
tin8_bkg0_ica.sampleinfo = tin8_bkg0.sampleinfo;
save([outdir 'tin8_bkg0_ica.mat'], 'tin8_bkg0_ica');

%tin8_bkg3
first = last+1;
last = first+length(tin8_bkg3.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin8_bkg3_ica = ft_selectdata(cfg, appended_ICA);
tin8_bkg3_ica.sampleinfo = tin8_bkg3.sampleinfo;
save([outdir 'tin8_bkg3_ica.mat'], 'tin8_bkg3_ica');

%tin8_bkg8
first = last+1;
last = first+length(tin8_bkg8.trial)-1;

cfg = [];
cfg.trials = (first:last);
tin8_bkg8_ica = ft_selectdata(cfg, appended_ICA);
tin8_bkg8_ica.sampleinfo = tin8_bkg8.sampleinfo;
save([outdir 'tin8_bkg8_ica.mat'], 'tin8_bkg8_ica');
