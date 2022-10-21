
%WIP
subinpath = '../mat_data/ICA/ID0905/';

for ii = 1:numel(conditions)

trig = cond.([conditions{ii} 'trig']);

triglabel = cond.([conditions{ii} 'label']);
nstim = numel(trig);

    for iii = 1:nstim

        %Only for pulse onset triggers (i.e GPP* PO* or PPP*)
        if any(regexp(triglabel{iii}, 'GPP_*')) || any(regexp(triglabel{iii}, 'PO_*')) || any(regexp(triglabel{iii}, 'PPP_*')) || any(regexp(triglabel{iii}, 'GO_*'))

            for i = 2%1:numel(sub_date.ID)
    
            fname = [(cond.(([conditions{ii} 'label'])){iii}) '_tlks'];
            fpath = ['../mat_data/timelockeds/ID' sub_date.ID{i} '/'];

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
    
            cfg = [];
            cfg.covariance = 'no';
            cfg.covariancewindow = 'prestim';
            cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"
            cfg.preproc.demean = 'yes';
            
            if any(regexp(triglabel{iii}, 'PO_*')) || any(regexp(triglabel{iii}, 'GO_*'))
                %Baseline window: 150ms before pulse/gap onset in PO/GO trials
                cfg.preproc.baselinewindow = [-0.150 0];
            elseif any(regexp(triglabel{iii}, 'PPP_*')) || any(regexp(triglabel{iii}, 'GPP_*'))
                %Baseline window: 150ms before gap/pp onset in GP/PP-trials
                cfg.preproc.baselinewindow = [-0.440 -0.290];
            end
            
            cfg.preproc.lpfilter = 'yes';
            cfg.preproc.lpfreq = 70;
            cfg.preproc.hpfilter = 'no';
            
            cfg.trials = tempdat.trialinfo == trig(iii);
            
            timelockeds = ft_timelockanalysis(cfg, tempdat);
            
            clear tempdat
            
            %save individual timelockeds
            save([fpath fname '.mat'], 'timelockeds');
                    
            cfg = [];
            timelockeds_cmb = ft_combineplanar(cfg, timelockeds);
    
            %save individual timelockeds with combined planar
            save([fpath fname '_cmb.mat'], 'timelockeds_cmb');
            
            clear timelockeds timelockeds_cmb;
            
            %For subjects
            end

        % if trigger of interest    
        end

    %For stim
    end

%For condition
end
