
%Frequancy analysis

%WIP
subinpath = '../mat_data/ICA/ID0905/';

for ii = 1:numel(conditions)

trig = cond.([conditions{ii} 'trig']);

triglabel = cond.([conditions{ii} 'label']);
nstim = numel(trig);

    for iii = 1:nstim

        %Only for pulse onset triggers (i.e GPP* PO* or PPP*)
        if any(regexp(triglabel{iii}, 'GPP_*')) || any(regexp(triglabel{iii}, 'PO_*')) %|| any(regexp(triglabel{iii}, 'PPP_*')) || any(regexp(triglabel{iii}, 'GO_*'))

            for i = 2%1:numel(sub_date.ID)
    
            fname = [(cond.(([conditions{ii} 'label'])){iii}) '_freq'];
            fpath = ['../mat_data/freqanalysis/ID' sub_date.ID{i} '/'];

            if ~exist(fpath, 'file');
            mkdir(fpath);
            end
            
            %check if file exist
            if exist([fpath fname '.mat'], 'file')
                warning([fname ' for subject: ID' sub_date.ID{i} ' exist'])
            continue
            end

            tempdat = load([subinpath conditions{ii} '_ica.mat']);
            tempdat = tempdat.([conditions{ii} '_ica']);
    
            cfg              = [];
            cfg.output       = 'powandcsd';
            cfg.channel      = 'MEG';
            cfg.method       = 'mtmconvol';
            cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
            cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
            cfg.toi          = -0.5:0.05:0.5;                  % time window "slides" from -0.5 to 0.5 sec in steps of 0.05 sec (50 ms)
            cfg.tapsmofrq    = 2;
            cfg.trials = tempdat.trialinfo == trig(iii);
            
            freqs = ft_freqanalysis(cfg, tempdat);
            
            %clear tempdat
            
            %save individual freqanalysis results
            save([fpath fname '.mat'], 'freqs');
                    
            %For subjects
            end

        % if trigger of interest    
        end

    %For stim
    end

%For condition
end

%% MultiplotTFR

freqs = load('../mat_data/freqanalysis/ID0905/PO_00_freq.mat');
freqs = freqs.freqs;

cfg = [];
freqs_cmb = ft_combineplanar(cfg, freqs);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.baseline = [-0.5 0];

ft_multiplotTFR(cfg, freqs_cmb);



