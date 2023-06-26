%% Process EOG

%Structure for subject timleockeds.avg
tlk_sub_eog = struct();

%For each condition
for ii = 1:length(temp_cond)

%Get triggers and n of triggers
trig = temp_stim.([temp_cond{ii} 'trig']);
nstim = numel(trig);

    for iii = 1:nstim

        subinpath = ['../processed_data/ICA/' 'ID' sub_date.ID{i} '/'];
        
        tempdat = load([subinpath temp_cond{ii} '_ica.mat']);
        tempdat = tempdat.([temp_cond{ii} '_ica']);

        cfg = [];
        cfg.covariance = 'no';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'yes'; % Keep all trials to clean later
        cfg.preproc.demean = 'yes';
        
        %Adjust baseline for experiment version and condition
        if exp_ver == "tinmeg1"
            if ismember(temp_cond{ii}, {'GO60', 'GO70', 'PO60', 'PO70'});
            %Baseline window: 150ms before pulse onset in PO/GO trials
            cfg.preproc.baselinewindow = [-0.150 0];
            elseif ismember(temp_cond{ii}, {'GP60', 'GP70'});
                %Baseline window variable timepoint: 150ms before gap onset in GP trials
                if trig(iii) == cond.tinmeg1.GP60trig(1) | trig(iii) == cond.tinmeg1.GP70trig(1) %ISI 0 GP60 | GP70
                cfg.preproc.baselinewindow = [-0.200 -0.050];
                elseif trig(iii) == cond.tinmeg1.GP60trig(2) | trig(iii) == cond.tinmeg1.GP70trig(2) %ISI 60
                    cfg.preproc.baselinewindow = [-0.260 -0.110];
                elseif trig(iii) == cond.tinmeg1.GP60trig(3) | trig(iii) == cond.tinmeg1.GP70trig(3) %ISI 120
                    cfg.preproc.baselinewindow = [-0.320 -0.170];
                elseif trig(iii) == cond.tinmeg1.GP60trig(4) | trig(iii) == cond.tinmeg1.GP70trig(4) %ISI 240
                    cfg.preproc.baselinewindow = [-0.440 -0.290];
                end
            end
        elseif exp_ver == "tinmeg2" || exp_ver == "tinmeg3"
            %Baseline window: 150ms before pulse onset
            cfg.preproc.baselinewindow = [-0.150 0];
            
            %Baseline window if preceded by Gap or pre-pulse
            if any(ismember(cond.tinmeg2.allGPP, trig(iii)))
                cfg.preproc.baselinewindow = [-0.440 -0.290];
            end

        end
        
        cfg.preproc.lpfilter = 'yes';
        cfg.preproc.lpfreq = 70;
        cfg.preproc.hpfilter = 'no';
        cfg.channel = {'EOG002'};
        cfg.trials = tempdat.trialinfo == trig(iii);
        
        tlk_sub_eog = ft_timelockanalysis(cfg, tempdat);
        
        clear tempdat
        
        %save individual timelockeds
        save([outdir cond.(sub_date.Exp{i}).([temp_cond{ii} 'label']){iii} '_EOG.mat'], 'tlk_sub_eog'); %clear tlk_sub
        
    end

end

clear ii iii trig nstim 

