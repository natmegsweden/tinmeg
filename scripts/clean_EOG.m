%% Clean EOG data

%check if subject is included
if ~any(ismember(EOG_clean.(temp_exp).ID, sub_date.ID{i}))

    %Save subject ID to first empty cell in struct
    EOG_clean.(temp_exp).ID{empty_IDx} = sub_date.ID{i};

    %For each condition
    for ii = 1:numel(cond.(sub_date.Exp{i}).all)
    
        %For each stim in condition
        for iii = 1:numel(cond.(sub_date.Exp{i}).([cond.(sub_date.Exp{i}).all{ii} 'label']))

            %Get condition and stimuli
            temp_cond = ([cond.(sub_date.Exp{i}).all{ii}]);
            temp_stim = cond.(sub_date.Exp{i}).([cond.(sub_date.Exp{i}).all{ii} 'label']){iii};

            %Write info
            disp(['Now showing EOG for ID' sub_date.ID{i} ': ' temp_cond '  |  ' temp_stim])
            
            %Get data for stim
            y = EOG_all_sub.(temp_exp).(temp_cond).(temp_stim){i}';

            %Get time
            t = (1:height(y))';
            
            %n trials
            n_trials = length(y);

            %name trials
            names = regexp(cellstr(sprintf('trial_n_%d ',1:n_trials)),' ','split');
            names = names{1,1};

            % Plot data
            figure('Position',[800 800 1000 800])
            hLines = plot(t, y, 'k');
            xline([101 101]);   %%%%%%%%%%%%%%%% not static plz
            xlim([1 201]);
            xticks([1:20:201]);
            xticklabels([-500:100:500]);
            
            % Start brushing mode and wait for user to hit "Enter" when done
            brush on
            disp('Hit Enter in comand window when done brushing')
            pause
            
            % Loop through each graphics object
            for k = 1:numel(hLines)
                % Check that the property is valid for that type of object
                % Also check if any points in that object are selected
                if isprop(hLines(k),'BrushData') && any(hLines(k).BrushData)
                    % Output the selected data to the base workspace with assigned name
                    ptsSelected = logical(hLines(k).BrushData.');
                    data = [t(ptsSelected) y(ptsSelected, k)];
                    assignin('base',names{k},data) %assign brushed data as variables names from "names"
                end
            end
            
            %Hacky way to get index of brushed variables
            vars = who('trial_n_*');
            clear(names{:});
            
            close;
            clear('colnum');

            for j = 1:numel(vars)
                colnum(j) = find(ismember(names, vars{j}));
            end
            
            %if any brushed trials, remove from temp variable y
            if exist('colnum','var') == 1 %if any removed, plot them as red
                y(:,[colnum]) = [];
            end
    
            % Write y to new struct
            EOG_clean.(temp_exp).(temp_cond).(temp_stim){empty_IDx} = y;

            clear names temp_cond temp_stim n_trials y
                
        %For stim
        end

    %For condition
    end

    clear temp_exp empty_IDx

%if subject already in struct
else
    disp([sub_date.ID{i} ': ' temp_cond '  |  ' temp_stim ' is already in struct'])
end
