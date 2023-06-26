%% Gather EOG

%Data path
inpath = ['../processed_data/timelockeds/' 'ID' sub_date.ID{i} '/EOG/'];

%Save subject ID to first empty cell in struct
EOG_all_sub.(temp_exp).ID{empty_IDx} = sub_date.ID{i};

%For each condition
for ii = 1:numel(cond.(sub_date.Exp{i}).all)

    %For each stim in condition
    for iii = 1:numel(cond.(sub_date.Exp{i}).([cond.(sub_date.Exp{i}).all{ii} 'label']))

        %Get condition and stimuli
        temp_cond = ([cond.(sub_date.Exp{i}).all{ii}]);
        temp_stim = cond.(sub_date.Exp{i}).([cond.(sub_date.Exp{i}).all{ii} 'label']){iii};

        %Specify filenames
        fname = [temp_stim '_EOG.mat'];

        disp(['Now gathering EOG for ID' sub_date.ID{i} ': ' temp_cond '  |  ' temp_stim])
        
        %load data
        tlk_sub_eog = load([inpath fname]);
        tlk_sub_eog = tlk_sub_eog.tlk_sub_eog;

        %Write time and sensor labels - NB! overwrites
        EOG_all_sub.(temp_exp).time = tlk_sub_eog.time;
        EOG_all_sub.(temp_exp).label = tlk_sub_eog.label;

        %Write average to first empty cell in structure
        EOG_all_sub.(temp_exp).(temp_cond).(temp_stim){empty_IDx} = tlk_sub_eog.trial(:,:);

        clear temp_cond temp_stim tlk_sub_eog
            
    %For stim
    end

%For condition
end

clear temp_exp empty_IDx