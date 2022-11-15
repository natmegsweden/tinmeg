
%Frequancy analysis

%WIP
subinpath = '../mat_data/ICA/ID0905/';

for ii = 1:3 %:numel(conditions)

trig = cond.([conditions{ii} 'trig']);

triglabel = cond.([conditions{ii} 'label']);
nstim = numel(trig);

    for iii = 1:nstim

        %Only for pulse onset triggers (i.e GPP* PO* or PPP*)
        if any(regexp(triglabel{iii}, 'GPP_*')) || any(regexp(triglabel{iii}, 'PO_*')) %|| any(regexp(triglabel{iii}, 'PPP_*')) || any(regexp(triglabel{iii}, 'GO_*'))
            
            %For subject
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
            cfg.foi          = [1:1:20];
            cfg.t_ftimwin    = repmat(0.1,1,20);
            cfg.toi          = '50%';
            cfg.tapsmofrq    = repmat(5,1,20);
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

freqsPO = load('../mat_data/freqanalysis/ID0905/PO_00_freq.mat');
freqsPO = freqsPO.freqs;

freqsGP = load('../mat_data/freqanalysis/ID0905/GPP_00_freq.mat');
freqsGP = freqsGP.freqs;

cfg = [];
freqsPO_cmb = ft_combineplanar(cfg, freqsPO);

cfg = [];
freqsGP_cmb = ft_combineplanar(cfg, freqsGP);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.baseline = [0 0.2];
cfg.title = 'PO';
ft_multiplotTFR(cfg, freqsPO_cmb);

cfg.title = 'GP';
cfg.baseline = [0 0.2];
ft_multiplotTFR(cfg, freqsGP_cmb);


cfg.channel = 'MEG1612+1613'
ft_singleplotTFR(cfg, freqsPO_cmb)
ax = gca;
ax.XLim = [0 0.2];
ax.Colormap = viridis(8)



freqsPO = load('../mat_data/freqanalysis/ID0905/PO_00_freq.mat');
freqsPO = freqsPO.freqs;

cfg = [];
freqsPO_cmb = ft_combineplanar(cfg, freqsPO);

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.baseline = [-0.4 0];
cfg.title = 'PO';
ft_multiplotTFR(cfg, freqsPO_cmb);



