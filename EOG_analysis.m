
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

%save(['../mat_data/epochs_eog.mat'], 'epochs_eog');
%load('../mat_data/epochs_eog.mat');


epochs_eog_resp = struct;

%Calculate max response for subjects
%NB, min between 75-150ms (sample x-x)
for i = 1:length(sub_date.ID)
    
    subinpath = [inpath 'ID' sub_date.ID{i} '/'];
    
    epochs_eog_resp.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    label = eval(['cond.' char(conditions(ii)) 'label']);
    
        for stim_index = 1:nstim
        
        %window:
        minwin = cell2mat(epochs_eog.(conditions{ii})(i,stim_index));
        %minresponse
        minresp = min(minwin(35:50));
        
        epochs_eog_resp.(conditions{ii}){i, stim_index} = minresp
            
        %For stim
        end
    
    %For conditions
    end

%For subjects
end


%To rectify vector a:
% a = [1 2 3 -4]
% a(a<0) = a(a<0)*-1

%[M, I] = min(mean(cell2mat(epochs_eog.PO60(:,6))));

%% PO60 plots

%Plot of mean EOG002 amplitude for all subjects PO60-stim
%NB "times" (X axis) defined manually, from C_process_timelockeds.m
figure
hold on
plot(mean(cell2mat(epochs_eog.PO60(:,1))));
plot(mean(cell2mat(epochs_eog.PO60(:,2))));
plot(mean(cell2mat(epochs_eog.PO60(:,3))));
plot(mean(cell2mat(epochs_eog.PO60(:,4))));
plot(mean(cell2mat(epochs_eog.PO60(:,5))));
plot(mean(cell2mat(epochs_eog.PO60(:,6))));

ax = gca;

plot([35 35], [ax.YLim(1), ax.YLim(2)], '--r');
plot([50 50], [ax.YLim(1), ax.YLim(2)], '--r');

legend('70', '75', '80', '85', '90', '95');
title('Average EOG amplitude after pulse only (n = 5)');

times = [-20:1:60]/200*1000;

ax = gca;
ax.XTickLabel = times(1:10:end);
ax.XLim = [0 80];

clear 'times';

xlabel('Latency (ms)');
ylabel('Amplitude (uV)');

[hleg,att] = legend('show');
title(hleg, 'Pulse level (dBA)');
hleg.Location = 'southwest';

hold off

% Plot average GP60/ISI
%NB "times" (X axis) defined manually, from C_process_timelockeds.m
figure
hold on
plot(mean(cell2mat(epochs_eog.GP60(:,1))));
plot(mean(cell2mat(epochs_eog.GP60(:,2))));
plot(mean(cell2mat(epochs_eog.GP60(:,3))));
plot(mean(cell2mat(epochs_eog.GP60(:,4))));

%plot(mean(cell2mat(epochs_eog.GO(:,1))), '--r', 'LineWidth', 1);

ax = gca;

plot([35 35], [ax.YLim(1), ax.YLim(2)], '--r');
plot([50 50], [ax.YLim(1), ax.YLim(2)], '--r');

times = [-20:1:60]/200*1000;

ax = gca;
ax.XTickLabel = times(1:10:end);
ax.XLim = [0 80];

legend('0','60','120','240');
title({'Average EOG amplitude with different', 'inter-stimulus interval (ISI) (n = 5)'});

[hleg,att] = legend('show');
title(hleg, 'ISI (ms)');
hleg.Location = 'southwest';

hold off

%Boxplot of maximum response amplitude PO60
boxlabel = [repmat({'70'}, length(epochs_eog_resp.PO60(:,1)), 1),
    repmat({'75'}, length(epochs_eog_resp.PO60(:,2)), 1),
    repmat({'80'}, length(epochs_eog_resp.PO60(:,3)), 1),
    repmat({'85'}, length(epochs_eog_resp.PO60(:,4)), 1),
    repmat({'90'}, length(epochs_eog_resp.PO60(:,5)), 1),
    repmat({'95'}, length(epochs_eog_resp.PO60(:,6)), 1)];

figure
hold on
boxplot([cell2mat(epochs_eog_resp.PO60(:,1)), 
    cell2mat(epochs_eog_resp.PO60(:,2)), 
    cell2mat(epochs_eog_resp.PO60(:,3)), 
    cell2mat(epochs_eog_resp.PO60(:,4)), 
    cell2mat(epochs_eog_resp.PO60(:,5)), 
    cell2mat(epochs_eog_resp.PO60(:,6))], boxlabel);

title('EOG amplitude after pulse only (n = 5)');

xlabel('Pulse level (dBA)');
ylabel('Amplitude (uV)');

hold off


%Boxplot Boxplot of maximum response amplitude GP60
boxlabel = [repmat({'ISI 0'}, length(epochs_eog_resp.GP60(:,1)), 1),
    repmat({'ISI 60'}, length(epochs_eog_resp.GP60(:,2)), 1),
    repmat({'ISI 120'}, length(epochs_eog_resp.GP60(:,3)), 1),
    repmat({'ISI 240'}, length(epochs_eog_resp.GP60(:,4)), 1)];

figure
hold on
boxplot([cell2mat(epochs_eog_resp.GP60(:,1)), 
    cell2mat(epochs_eog_resp.GP60(:,2)), 
    cell2mat(epochs_eog_resp.GP60(:,3)), 
    cell2mat(epochs_eog_resp.GP60(:,4))], boxlabel);

title('EOG amplitude after gap-pulse (n = 5)');

xlabel('Inter-stimulus interval (ms)');
ylabel('Amplitude (uV)');

hold off
