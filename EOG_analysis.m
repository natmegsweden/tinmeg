
%TO DO, limit time-window for min() (i.e max EOG-response)

epochs_eog = struct;

inpath = '../mat_data/timelockeds/';

%load EOG data to structure
for i = 1:length(sub_date.ID)
    
    subinpath = [inpath 'ID' sub_date.ID{i} '/'];
    
    epochs_eog.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    label = eval(['cond.' char(conditions(ii)) 'label']);
    
        for stim_index = 1:nstim    

        VAR = 'eog_timelockeds';
        T = load([subinpath char(label(stim_index)) '_eog' '.mat'], VAR);
        T = T.(VAR).avg(2,:);

        epochs_eog.(conditions{ii}){i, stim_index} = T
        
        clear('T', 'VAR')
        
        %For stim
        end
    
    %For conditions
    end

%For subjects
end

epochs_eog_resp = struct;

%Calculate max response for subjects
for i = 1:length(sub_date.ID)
    
    subinpath = [inpath 'ID' sub_date.ID{i} '/'];
    
    epochs_eog_resp.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    label = eval(['cond.' char(conditions(ii)) 'label']);
    
        for stim_index = 1:nstim
        
        %TO DO: define time-window
        minresp = min(cell2mat(epochs_eog.(conditions{ii})(i,stim_index)));
        
        epochs_eog_resp.(conditions{ii}){i, stim_index} = minresp
            
        %For stim
        end
    
    %For conditions
    end

%For subjects
end

%save(['../mat_data/epochs_eog.mat'], 'epochs_eog');
%load('../mat_data/epochs_eog.mat');

figure
hold on
boxplot([epochs_eog_resp.PO60{1,:},
    epochs_eog_resp.PO60{2,:},
    epochs_eog_resp.PO60{3,:},
    epochs_eog_resp.PO60{4,:},
    epochs_eog_resp.PO60{5,:},
    epochs_eog_resp.PO60{6,:}], 'Notch', 'on', 'Labels', {'70', '75', '80', '85'})

hold off

%Boxplot of maximum response amplitude for subjects 
boxlabel = [repmat({'70'}, length(epochs_eog_resp.PO60(:,1)), 1),
    repmat({'75'}, length(epochs_eog_resp.PO60(:,2)), 1),
    repmat({'80'}, length(epochs_eog_resp.PO60(:,3)), 1),
    repmat({'85'}, length(epochs_eog_resp.PO60(:,4)), 1),
    repmat({'90'}, length(epochs_eog_resp.PO60(:,5)), 1),
    repmat({'95'}, length(epochs_eog_resp.PO60(:,6)), 1)]

figure
hold on
boxplot([cell2mat(epochs_eog_resp.PO60(:,1)), 
    cell2mat(epochs_eog_resp.PO60(:,2)), 
    cell2mat(epochs_eog_resp.PO60(:,3)), 
    cell2mat(epochs_eog_resp.PO60(:,4)), 
    cell2mat(epochs_eog_resp.PO60(:,5)), 
    cell2mat(epochs_eog_resp.PO60(:,6))], boxlabel)

title('EOG amplitude after pulse only (n = 5)');

xlabel('Pulse level');
ylabel('Amplitude (uV)');

hold off


%Plot of mean EOG002 amplitude for all subjects PO60-stim
figure
hold on
plot(mean(cell2mat(epochs_eog.PO60(:,1))));
plot(mean(cell2mat(epochs_eog.PO60(:,2))));
plot(mean(cell2mat(epochs_eog.PO60(:,3))));
plot(mean(cell2mat(epochs_eog.PO60(:,4))));
plot(mean(cell2mat(epochs_eog.PO60(:,5))));
plot(mean(cell2mat(epochs_eog.PO60(:,6))));

legend('70', '75', '80', '85', '90', '95');
title('Average EOG amplitude after pulse only (n = 5)');

xlabel('Latency (2ms sample)');
ylabel('Amplitude (uV)');

[hleg,att] = legend('show');
title(hleg, 'Pulse level (dBA)');

hold off

%To rectify vector a:
% a = [1 2 3 -4]
% a(a<0) = a(a<0)*-1

[M, I] = min(mean(cell2mat(epochs_eog.PO60(:,6))));

