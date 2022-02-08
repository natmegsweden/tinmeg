
% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

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

%% ICA


%% Load RS data and create pseudo epochs

for i = 1%:4%length(sub_date.ID);

    outdir = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/'];
    
    fname = ['rawRS_meg_dsc.mat'];
    fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

    %check if file exist
    if exist(fpath, 'file')
    warning([fname ' for subject: ID' sub_date.ID{i} ' exist'])
    continue
    end
    
    %load
    rawRS_meg_ds = load(['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/rawRS_meg_ds.mat']);
    rawRS_meg_ds = rawRS_meg_ds.rawRS_meg_ds;
    
    cfg = [];
    cfg.channel = 'meg';
    data = ft_selectdata(cfg, rawRS_meg_ds);
    
    %Make pseudo-epochs
    cfg = [];
    cfg.length  = 2;
    cfg.overlap = 0.5;       % 50% overlap
    epo = ft_redefinetrial(cfg, data);

    %https://www.fieldtriptoolbox.org/faq/how_can_i_do_time-frequency_analysis_on_continuous_data/
    cfg = [];
    cfg.method     = 'mtmfft'
    %cfg.taper      = 'hanning'
    cfg.tapsmofrq  = 2;
    cfg.foilim     = [1 30];
    cfg.keeptrials = 'no';
    cfg.output     = 'fooof_aperiodic';
    freq_segmented = ft_freqanalysis(cfg, epo)
    
    begsample = epo.sampleinfo(:,1);
    endsample = epo.sampleinfo(:,2);
    time = ((begsample+endsample)/2) / epo.fsample;

    freq_continuous           = freq_segmented;
    freq_continuous.powspctrm = permute(freq_segmented.powspctrm, [2, 3, 1]);
    freq_continuous.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
    freq_continuous.time      = time;             % add the description of the time dimension

    %https://www.fieldtriptoolbox.org/faq/how_can_i_do_time-frequency_analysis_on_continuous_data/
    
    plot(freq_segmented.freq, mean(freq_segmented.powspctrm, 1));
    
    
    






    
    %ft_rejectvisual
    cfg = [];
    cfg.method = 'summary';
    cfg.keepchannel = 'yes';
    cfg.channel = 'MEGMAG';
    cfg.layout = 'neuromag306all.lay';

    rawRS_meg_dsc = ft_rejectvisual(cfg,epo);

    cfg.channel = 'MEGGRAD';

    rawRS_meg_dsc = ft_rejectvisual(cfg, rawRS_meg_dsc);
    
    save([outdir 'rawRS_meg_dsc'], 'rawRS_meg_dsc')
    
    %Load if FOOOF does not exist?

    
    %Calculate FOOOF (https://www.fieldtriptoolbox.org/example/fooof/)
    cfg               = [];
    cfg.foilim        = [1 100];
    %cfg.pad           = 4;
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    %cfg.output        = 'fooof_aperiodic';
    %fractal = ft_freqanalysis(cfg, epo);
    cfg.output        = 'pow';
    original = ft_freqanalysis(cfg, epo);

end

cfg = [];

ft_singleplotTFR(cfg, original);

