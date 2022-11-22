
% Run setup for common variables if not already loaded
if exist('sub_date', 'var') == 0; run Conditions_triggers.m; end;

%% Process timelockeds per condition, gather subject - and grand average

%Structure for subject timleockeds.avg
tlk_sub = struct();

%Structure for subject timleockeds.avg with combined planar gradiometers
tlk_sub_cmb = struct();

for ii = 1:length(conditions)

trig = cond.([conditions{ii} 'trig']);
nstim = numel(trig);

    for iii = 1:nstim

        for i = 1:numel(sub_date.ID)

        subinpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
        
        tempdat = load([subinpath conditions{ii} 'ica.mat']);
        tempdat = tempdat.([conditions{ii} 'ica']);

        cfg = [];
        cfg.covariance = 'no';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"
        cfg.preproc.demean = 'yes';
        
        if ismember(conditions{ii}, {'GO60', 'GO70', 'PO60', 'PO70'});
        %Baseline window: 150ms before pulse onset in PO trials
        cfg.preproc.baselinewindow = [-0.150 0];
        elseif ismember(conditions{ii}, {'GP60', 'GP70'});
            %Baseline window variable timepoint: 150ms before gap onset in GP trials
            if iii == 1 %ISI 0
            cfg.preproc.baselinewindow = [-0.200 -0.050];
            elseif iii == 2 %ISI 60
                cfg.preproc.baselinewindow = [-0.260 -0.110];
            elseif iii == 3 %ISI 120
                cfg.preproc.baselinewindow = [-0.320 -0.170];
            elseif iii == 4 %ISI 240
                cfg.preproc.baselinewindow = [-0.440 -0.290];
            end
        end
        
        cfg.preproc.lpfilter = 'yes';
        cfg.preproc.lpfreq = 70;
        cfg.preproc.hpfilter = 'no';
        cfg.channel = 'MEG';
        cfg.trials = tempdat.trialinfo == trig(iii);
        
        timelockeds = ft_timelockanalysis(cfg, tempdat);
        
        clear tempdat
        
        %save individual timelockeds
        save(['../mat_data/timelockeds/ID' sub_date.ID{i} '/' (cond.(([conditions{ii} 'label'])){iii}) '_tlks.mat'], 'timelockeds');
                
        cfg = [];
        timelockeds_cmb = ft_combineplanar(cfg, timelockeds);

        %save individual timelockeds with combined planar
        save(['../mat_data/timelockeds/ID' sub_date.ID{i} '/' (cond.(([conditions{ii} 'label'])){iii}) '_tlks_cmb.mat'], 'timelockeds_cmb');
        
        %Write averages to structs for easy access
        tlk_sub_cmb.(conditions{ii}){i, iii} = timelockeds_cmb.avg;
        tlk_sub.(conditions{ii}){i, iii} = timelockeds.avg;
        
        end
        
    end

end

clear i ii iii trig nstim

save('../mat_data/timelockeds/tlk_sub.mat', 'tlk_sub');
save('../mat_data/timelockeds/tlk_sub_cmb.mat', 'tlk_sub_cmb');

