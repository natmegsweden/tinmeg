%% Run ICA on downsampled data to remove heartbeat artefacts. NB: NOT EOG - high correlation with response excpected

%Check if subject dir exist, skip or create folder
if ~exist(outdir, 'file');
    mkdir(outdir);
end

%Check if figure dir for sub exist, create/define
if ~exist(figoutdir, 'file');
    mkdir(figoutdir);
end

%Find cleaned up files (different filenames for TinMEG 1/2/3
if sub_date.Exp{i} == 'tinmeg1'
    subpath = [inpath '/'];
    fnames = find_files(subpath, {'_cl'}, 'ds');
elseif sub_date.Exp{i} == 'tinmeg2'
    subpath = [inpath '/'];
    fnames = find_files(subpath, {'_cl'}, 'ds');
elseif sub_date.Exp{i} == 'tinmeg3'
    subpath = [inpath '/'];
    fnames = find_files(subpath, {'_cl'}, 'ds');
end

all_dat = cell(1);

%Load each file
for iii = 1:numel(fnames)

    all_dat{iii} = importdata([inpath fnames{iii}]);

end

%Append data
cfg = [];
cfg.keepsampleinfo = 'no';
appended = ft_appenddata(cfg, all_dat{1:end});

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

saveas(gcf, [figoutdir 'comp_figure.pdf'])

close(1)

% Remove components
cfg = [];
cfg.component   = ecg_comp_idx; %index of components to remove
cfg.channel     = 'MEG';
cfg.updatesens  = 'no'; %
appended_ICA = ft_rejectcomponent(cfg, comp, appended);


%Split conditions and save

%NB, order from fnames (alphabetical):
%GO, GP60, GP70, PO60, PO70

if sub_date.Exp{i} == 'tinmeg1'
    
    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg1.PO60trig);
    PO60_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'PO60_ica.mat'], 'PO60_ica'); clear PO60_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg1.PO70trig);
    PO70_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'PO70_ica.mat'], 'PO70_ica'); clear PO70_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg1.GP60trig);
    GP60_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'GP60_ica.mat'], 'GP60_ica'); clear GP60_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg1.GP70trig);
    GP70_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'GP70_ica.mat'], 'GP70_ica'); clear GP70_ica;
    
    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg1.GOtrig);
    GO_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'GO_ica.mat'], 'GO_ica'); clear GO_ica;

elseif sub_date.Exp{i} == 'tinmeg2'
%Tin0
    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin0_bkg0trig);
    tin0_bkg0_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin0_bkg0_ica.mat'], 'tin0_bkg0_ica'); clear tin0_bkg0_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin0_bkg3trig);
    tin0_bkg3_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin0_bkg3_ica.mat'], 'tin0_bkg3_ica'); clear tin0_bkg3_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin0_bkg8trig);
    tin0_bkg8_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin0_bkg8_ica.mat'], 'tin0_bkg8_ica'); clear tin0_bkg8_ica;
%Tin3
    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin3_bkg0trig);
    tin3_bkg0_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin3_bkg0_ica.mat'], 'tin3_bkg0_ica'); clear tin3_bkg0_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin3_bkg3trig);
    tin3_bkg3_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin3_bkg3_ica.mat'], 'tin3_bkg3_ica'); clear tin3_bkg3_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin3_bkg8trig);
    tin3_bkg8_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin3_bkg8_ica.mat'], 'tin3_bkg8_ica'); clear tin3_bkg8_ica;
%Tin8
    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin8_bkg0trig);
    tin8_bkg0_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin8_bkg0_ica.mat'], 'tin8_bkg0_ica'); clear tin8_bkg0_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin8_bkg3trig);
    tin8_bkg3_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin8_bkg3_ica.mat'], 'tin8_bkg3_ica'); clear tin8_bkg3_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg2.tin8_bkg8trig);
    tin8_bkg8_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin8_bkg8_ica.mat'], 'tin8_bkg8_ica'); clear tin8_bkg8_ica;

elseif sub_date.Exp{i} == 'tinmeg3'

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg3.tin0_bkg0trig);
    tin0_bkg0_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin0_bkg0_ica.mat'], 'tin0_bkg0_ica'); clear tin0_bkg0_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg3.tin0_bkg3trig);
    tin0_bkg3_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin0_bkg3_ica.mat'], 'tin0_bkg3_ica'); clear tin0_bkg3_ica;

    cfg = [];
    cfg.trials = ismember(appended_ICA.trialinfo, cond.tinmeg3.tin0_bkg8trig);
    tin0_bkg8_ica = ft_selectdata(cfg, appended_ICA);
    save([outdir 'tin0_bkg8_ica.mat'], 'tin0_bkg8_ica'); clear tin0_bkg8_ica;

end

clear appended_ICA appended temp data_ecg comp_ecg freq fd_comp fdcomp comp ecg;