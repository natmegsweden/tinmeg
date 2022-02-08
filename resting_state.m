
%Local fieldtrip path
addpath D:\MATLAB\fieldtrip-20220206;
ft_defaults;

sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

%%

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%Sample rate to downsample to
fs_ds = 200;

% Create cell array for subjects filepaths
subpaths = cell(1);

for i = 1%:2%length(sub_date.ID);

    outdir = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/'];

    % Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
    subpath = [meg_data_path 'NatMEG_' sub_date.ID{i} '/' sub_date.date{i} '/'];
    fnames = find_files(subpath, {'restingstate', 'tsss'}, {'ds', 'test'});

    subpaths{i,1} = sub_date.ID{i}; %Include ID for tracking

        for fileindex = 1:length(fnames);
            subpaths{i,1+fileindex} = [subpath char(fnames(fileindex))]; % NB! n of files differ between rows, some subjects have empty columns
        end

    clear fnames subpath fileindex

    sub_data = subpaths(i,2:max(find(~cellfun(@isempty,subpaths(i,:)))));

    %Preprocess
    cfg = [];
    cfg.dataset     = sub_data;
    cfg.continuous  = 'yes';
    cfg.demean      = 'yes';

    rawRS_meg = ft_preprocessing(cfg);
    
    save([outdir 'rawRS_meg'], 'rawRS_meg')

    %downsample and save
    cfg = [];
    cfg.resamplefs = fs_ds;

    rawRS_meg_ds = ft_resampledata(cfg, rawRS_meg);

    save([outdir 'rawRS_meg_ds'], 'rawRS_meg_ds')
    
end

%% Load RS data, clean and create pseudo epochs

for i = 1:length(sub_date.ID);

    outdir = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/'];
    
    fname = ['rawRS_meg_dsc.mat'];
    fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

    %check if file exist
    if exist(fpath, 'file')
    warning([fname ' for subject: ID' sub_date.ID{i} ' exist - skipping reject_visual']);
    
    rawRS_meg_dsc = load(fpath);
    rawRS_meg_dsc = rawRS_meg_dsc.rawRS_meg_dsc;
   
    elseif ~exist(fpath, 'file')
        
    %load
    rawRS_meg_ds = load(['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/rawRS_meg_ds.mat']);
    rawRS_meg_ds = rawRS_meg_ds.rawRS_meg_ds;
    
    %Make pseudo-epochs
    cfg = [];
    cfg.length  = 2;
    cfg.overlap = 0.5;       % 50% overlap
    epo = ft_redefinetrial(cfg, rawRS_meg_ds);
    
    %ft_rejectvisual
    cfg = [];
    cfg.method = 'summary';
    cfg.keepchannel = 'yes';
    cfg.channel = 'MEGMAG';
    cfg.layout = 'neuromag306all.lay';

    rawRS_meg_dsc = ft_rejectvisual(cfg,epo);
    
    clear epo
    
    cfg.channel = 'MEGGRAD';

    rawRS_meg_dsc = ft_rejectvisual(cfg, rawRS_meg_dsc);
    
    save(fpath, 'rawRS_meg_dsc');

    end

    cfg = [];
    cfg.channel = 'meg';
    meg_data = ft_selectdata(cfg, rawRS_meg_dsc);

    % compute the fractal and original spectra
    cfg               = [];
    cfg.foilim        = [2 30];
    cfg.taper         = 'hanning';
    cfg.pad           = 'nextpow2';
    cfg.tapsmofrq     = 0;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_aperiodic';
    fractal = ft_freqanalysis(cfg, meg_data);
    cfg.output        = 'pow';
    original = ft_freqanalysis(cfg, meg_data);
    
    clear meg_data
    
    % subtract the fractal component from the power spectrum
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2-x1';
    oscillatory = ft_math(cfg, fractal, original);
    
    figure; hold on;
    title(sub_date.ID{i});
    plot(fractal.freq, mean(fractal.powspctrm,1));
    plot(original.freq, mean(original.powspctrm,1));
    plot(oscillatory.freq, mean(oscillatory.powspctrm,1));
    saveas(gcf, ['../output/' sub_date.ID{i} '_fooof.png']);
    close
    
    original_pwspc{i} = mean(original.powspctrm,1);
    fract_pwspc{i} = mean(fractal.powspctrm,1);
    osc_pwspc{i} = mean(oscillatory.powspctrm, 1);
    
    freqs = oscillatory.freq;

end

% save('../output/osc_pwspc.mat', 'osc_pwspc');
% save('../output/original_pwspc.mat', 'original_pwspc');
% save('../output/osc_fract.mat', 'fract_pwspc');

%%
for i = 1:22
   
    [M, I] = max(osc_pwspc{i}(16:27))
    
    maxPSDfreq(i) = freqs(16+I-1);
    
end

%% Plots

load('../output/osc_pwspc.mat');
load('../output/original_pwspc.mat');
load('../output/osc_fract.mat');

blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];

figure('Units', 'centimeters', 'Position', [5 5 15 10]);
histogram(maxPSDfreq)
set(gca,'FontSize',10)
%title('Distribution of maximum PSD frequency in the Alpha band');
xlabel('Maximum PSD frequency (Hz)');
ylabel('n count');
ylim([0 10]);

saveas(gcf, ['../output/RS_histo.pdf']);

%FoooF demo
figure('Units', 'centimeters', 'Position', [5 5 15 12]);
subplot(2,2,1); hold on;
set(gca,'FontSize',10)
plot(freqs, osc_pwspc{2});
plot(freqs, fract_pwspc{2});
plot(freqs, original_pwspc{2});
ylabel('PSD (AU)');
ylim([-1*10^-26 20*10^-26]);
xlim([1 30]);

subplot(2,2,2); hold on;
set(gca,'FontSize',10)
plot(freqs, osc_pwspc{5});
plot(freqs, fract_pwspc{5});
plot(freqs, original_pwspc{5});
legend('Periodic', 'Aperiodic', 'Original');
ylim([-1*10^-26 10*10^-25]);
xlim([1 30]);

subplot(2,2,3); hold on;
set(gca,'FontSize',10)
plot(freqs, osc_pwspc{18});
plot(freqs, fract_pwspc{18});
plot(freqs, original_pwspc{18});
ylabel('PSD (AU)');
xlabel('Frequency (Hz)');
ylim([-1*10^-26 20*10^-26]);
xlim([1 30]);

subplot(2,2,4); hold on;
set(gca,'FontSize',10)
plot(freqs, osc_pwspc{22});
plot(freqs, fract_pwspc{22});
plot(freqs, original_pwspc{22});
xlabel('Frequency (Hz)');
ylim([-1*10^-26 20*10^-26]);
xlim([1 30]);

saveas(gcf, ['../output/RS_subplots.pdf']);

%Summary plot
figure('Units', 'centimeters', 'Position', [5 5 15 8]); hold on;
set(gca,'FontSize',10)
for i = [1:10 12:18 20 21 22]
    if i < 22
    plot(freqs, osc_pwspc{i}, 'Color', [0 0 0 0.5], 'HandleVisibility','off')
    else
    plot(freqs, osc_pwspc{i}, 'Color', [0 0 0 0.5])
    end
    tempmean(i,1:73) = osc_pwspc{i};
end
plot(freqs, osc_pwspc{18}, 'Color', [0 0 0], 'LineWidth', 1.5)
plot(freqs, mean(tempmean), 'Color', [0.75 0 0], 'LineWidth', 1.5);
ylim([-1*10^-25 20*10^-25]);
xlim([1 20]);
xlabel('Frequency (Hz)');
ylabel('PSD (AU)');
legend('Subject 1 - 21', 'Subject 22', 'Average PSD');

patch('Faces', [1 2 3 4], 'Vertices', [3.5 -1*10^-25; 3.5 20*10^-25; 7 20*10^-25; 7 -1*10^-25], ...
'FaceColor', blue, 'FaceAlpha', 0.1, 'EdgeAlpha', 0, 'DisplayName', 'Theta (3.5 - 7 Hz)');
patch('Faces', [1 2 3 4], 'Vertices', [8 -1*10^-25; 8 20*10^-25; 12 20*10^-25; 12 -1*10^-25], ...
'FaceColor', orange, 'FaceAlpha', 0.1, 'EdgeAlpha', 0, 'DisplayName', 'Alpha (8 - 12 Hz)');

saveas(gcf, ['../output/RS_average.pdf']);