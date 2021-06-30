
%epochs_eog.(conditions{ii}){i, stim_index} = ft_timelockanalysis(cfg, temptrials);

epochs_eog = struct;

inpath = '../mat_data/timelockeds/'

for i = 1:length(sub_date.ID)
    
    subinpath = [inpath 'ID' sub_date.ID{i} '/'];
    
    epochs_eog.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    label = eval(['cond.' char(conditions(ii)) 'label']);
    
        for stim_index = 1:nstim    

        VAR = 'eog_timelockeds'
        T = load([subinpath char(label(stim_index)) '_eog' '.mat'], VAR);
        T = T.(VAR)

        epochs_eog.(conditions{ii}){i, stim_index} = T

        %epochs_eog.(conditions{ii}){i, stim_index} = ft_timelockanalysis(cfg, temptrials);
        end
        
    end
    
end