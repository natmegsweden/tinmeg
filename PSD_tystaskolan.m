%% Read single subject resting state and create whole brain PSD figure for tysta skolan grant application 2021

%Setup paths
meg_data_path = '/archive/20061_tinnitus/MEG/NatMEG_0539/210208/';
           
% List of all filenames to import                
filenames = {...
        'restingstate_mc_avgtrans_ds_5_tsss_corr98.fif' 
        'restingstate-1_mc_avgtrans_ds_5_tsss_corr98.fif'
            }; 
        
% Define where to put output data
output_path = '/home/nikedv/TinMEG1/Analysis Output';

%% Preprocess and append fif

data_meg = cell(1,(length(filenames)));

for fileindex = 1:length(filenames);
    
cfg = [];
cfg.dataset     = [meg_data_path filenames{fileindex}];
cfg.continuous  = 'yes';
cfg.demean      = 'yes';

data_meg{fileindex}  = ft_preprocessing(cfg);

end

cfg = [];
data_meg_append = ft_appenddata(cfg, data_meg{:});

%% Preprocess single fif

cfg = [];
cfg.dataset     = [meg_data_path filenames{1}];
cfg.continuous  = 'yes';
cfg.demean      = 'yes';

data_meg_single  = ft_preprocessing(cfg);

%% 

%Using data_meg_append gives error: Reference to non-existent field 'chanunit'.

%Make pseudo-epochs
cfg = [];
cfg.length  = 3;
cfg.overlap = .5;       % 50% overlap
epo = ft_redefinetrial(cfg, data_meg_single);       

% Get PSD
cfg = [];
cfg.channel = 'meg';          % Only MEG channels
cfg.method  = 'mtmfft';
cfg.output  = 'pow';
cfg.taper   = 'hanning';      % Single taper
cfg.foilim  = [1 48];
cfg.pad     = 'nextpow2';

pow = ft_freqanalysis(cfg, epo);

%Plot power spectrum
cfg = [];
cfg.linewidth = 2;
ft_singleplotER(cfg, pow)
%xlim([0 40])
ylim([0 15*10^-26])
title('Whole brain PSD')
