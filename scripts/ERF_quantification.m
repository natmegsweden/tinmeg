%% Quantify ERF responses from timelockeds

if exist('tlk_all_sub','var') == 0
load('../processed_data/timelockeds/aggregated/tlk_all_sub_cmb.mat', 'tlk_all_sub');
end

%Load names of highest N1 ERF gradiometers in tinmeg1 PO60_90
load('../analysis_output/vars/topgrads_tinemg1_PO60_90_N1.mat');

%Create empty structure
ERFs_n1 = struct();

%Get list of experiments
exps = unique(sub_date.Exp);

%Decide time window of interest
toi_start = 0.050; % start (seconds)
toi_end = 0.150;   % end (seconds)

%for each experiment
for i = 1:numel(exps)

    %Find index for TOI
    toi1 = find(tlk_all_sub.(exps{i}).time == toi_start);
    toi2 = find(tlk_all_sub.(exps{i}).time == toi_end);

    %Find index of top responding grad
    topLind = find(ismember(tlk_all_sub.(exps{i}).label,topL));
    topRind = find(ismember(tlk_all_sub.(exps{i}).label,topR));

    %For each experimental condition
    for ii = 1:numel(cond.(exps{i}).all)

        exp_conds = cond.(exps{i}).all;

        %For each stim in condition
        for iii = 1:numel(cond.(exps{i}).([cond.(exps{i}).all{ii} 'label']));

            temp_stims = cond.(exps{i}).([cond.(exps{i}).all{ii} 'label']);
    
            %For each subject in experiment
            for iiii = 1:numel(tlk_all_sub.(exps{i}).ID)

                %Left side
                ERFs_n1.(exps{i}).(exp_conds{ii}).L.(temp_stims{iii}){iiii} = mean(tlk_all_sub.(exps{i}).(exp_conds{ii}).(temp_stims{iii}){iiii}(topLind,toi1:toi2))

                %Right side
                ERFs_n1.(exps{i}).(exp_conds{ii}).R.(temp_stims{iii}){iiii} = mean(tlk_all_sub.(exps{i}).(exp_conds{ii}).(temp_stims{iii}){iiii}(topRind,toi1:toi2))
    
            end

        end

    end

end

save('../analysis_output/ERF/top_grad_N1.mat', 'ERFs_n1');