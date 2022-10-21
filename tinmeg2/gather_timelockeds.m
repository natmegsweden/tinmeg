%% SHOULD load if exist

%Structure for subject timleockeds.avg
tlk_sub = struct();

%Structure for subject timleockeds.avg with combined planar gradiometers
tlk_sub_cmb = struct();

%% Gather timelockeds

for ii = 1:numel(conditions)

trig = cond.([conditions{ii} 'trig']);
triglabel = cond.([conditions{ii} 'label']);
nstim = numel(trig);

    for iii = 1:nstim
            
            %Triggers of interest (i.e. GPG is redundandt)
            if any(regexp(triglabel{iii}, 'GPP_*')) || any(regexp(triglabel{iii}, 'PO_*')) || any(regexp(triglabel{iii}, 'PPP_*')) || any(regexp(triglabel{iii}, 'GO_*'))

            for i = 1:numel(sub_date.ID)

            fname = [(cond.(([conditions{ii} 'label'])){iii}) '_tlks'];
            fpath = ['../mat_data/timelockeds/ID' sub_date.ID{i} '/'];
            
            load([fpath fname '.mat']);
            load([fpath fname '_cmb.mat']);
            
            %Write averages to structs for easy access
            tlk_sub.(conditions{ii}){i, iii} = timelockeds.avg;
            tlk_sub.ID(i) = sub_date.ID{i};

            tlk_sub_cmb.(conditions{ii}){i, iii} = timelockeds_cmb.avg;
            tlk_sub_cmb.ID(i) = sub_date.ID{i};

            disp(['Gathered ' conditions{ii} ' for ID' num2str(sub_date.ID{i})]);

            %For subjects
            end
            
            %if trigger of interest
            end

    %For stim
    end

%For condition
end

save('../mat_data/timelockeds/tinmeg2/tlk_sub.mat', 'tlk_sub');
save('../mat_data/timelockeds/tinmeg2/tlk_sub_cmb.mat', 'tlk_sub_cmb');
