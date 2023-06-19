%% Gather timelockeds

%Create empty structure
tlk_all_sub = struct();

%For each subject
for i = 1:numel(sub_date.ID)

    %Data path
    inpath = ['../processed_data/timelockeds/' 'ID' sub_date.ID{i} '/'];
    
    %Get experiment version
    temp_exp = (sub_date.Exp{i});

    IDs = find(~cellfun(@isempty, tlk_all_sub.(temp_exp).ID));

    %Find first empty cell for ID
    empty_IDx = 1 + IDs(end);

    %check if subject is included
    if ~any(ismember(tlk_all_sub.(temp_exp).ID, sub_date.ID{i}))

        %Save subject ID to first empty cell in struct
        tlk_all_sub.(temp_exp).ID{empty_IDx} = sub_date.ID{i};
    
        %For each condition
        for ii = 1:numel(cond.(sub_date.Exp{i}).all)
        
            %For each stim in condition
            for iii = 1:numel(cond.(sub_date.Exp{i}).([cond.(sub_date.Exp{i}).all{ii} 'label']))
    
                %Get condition and stimuli
                temp_cond = ([cond.(sub_date.Exp{i}).all{ii}]);
                temp_stim = cond.(sub_date.Exp{i}).([cond.(sub_date.Exp{i}).all{ii} 'label']){iii};
    
                %Specify filenames
                fname = [temp_stim '_tlk.mat'];
                fname_cmb = [temp_stim '_tlk_cmb.mat'];
    
                disp(['Now gathering tlks for ID' sub_date.ID{i} ': ' temp_cond '  |  ' temp_stim])
                
                %load data
                tlk_sub_cmb = load([inpath fname_cmb]);
                tlk_sub_cmb = tlk_sub_cmb.tlk_sub_cmb;
    
                %Write time and sensor labels - NB! overwrites
                tlk_all_sub.(temp_exp).time = tlk_sub_cmb.time;
                tlk_all_sub.(temp_exp).label = tlk_sub_cmb.label;
    
                %Write average to first empty cell in structure
                tlk_all_sub.(temp_exp).(temp_cond).(temp_stim){empty_IDx} = tlk_sub_cmb.avg;
    
                clear temp_cond temp_stim tlk_sub tlk_sub_cmb
                    
            %For stim
            end
    
        %For condition
        end
    
        clear temp_exp empty_IDx

    %if subject already in struct
    else
        disp([sub_date.ID{i} ': ' temp_cond '  |  ' temp_stim ' is already in struct'])
    end

%For subject
end

save('../processed_data/timelockeds/aggregated/tlk_all_sub_cmb.mat', 'tlk_all_sub');
