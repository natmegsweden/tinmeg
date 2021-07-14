
inpath = '../mat_data/timelockeds/';

load('../mat_data/epochs_eog.mat');
load('../mat_data/epochs_eog_resp.mat');
load('../mat_data/epochs_eog_all.mat');
load(['../mat_data/epochs_eog_avgrast.mat');


%% load EOG data to structure

epochs_eog = struct;

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

%% load EOG keeptrials data to structure

epochs_eog_all = struct;

for i = 1:length(sub_date.ID)
    
    subinpath = [inpath 'ID' sub_date.ID{i} '/'];
    
    epochs_eog_all.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    label = eval(['cond.' char(conditions(ii)) 'label']);
    
        for stim_index = 1:nstim

        VAR = 'eog_timelockeds_all';
        T = load([subinpath char(label(stim_index)) '_eog_all' '.mat'], VAR);
        T = T.(VAR).trial(:,:);

        epochs_eog_all.(conditions{ii}){i, stim_index} = T
        
        clear('T', 'VAR')
        
        %For stim
        end
    
    %For conditions
    end

%For subjects
end
    
%save(['../mat_data/epochs_eog_all.mat'], 'epochs_eog_all');


%% find biggest response from structure

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

%save(['../mat_data/epochs_eog_resp.mat'], 'epochs_eog_resp');


%To rectify vector a:
% a = [1 2 3 -4]
% a(a<0) = a(a<0)*-1

%[M, I] = min(mean(cell2mat(epochs_eog.PO60(:,6))));

%% Create struct of average rasters

epochs_eog_avgrast = struct;

subinpath = [inpath 'ID' sub_date.ID{i} '/'];
epochs_eog_avgrast.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    label = eval(['cond.' char(conditions(ii)) 'label']);
    
        for stim_index = 1:nstim
                
                %sum up matrices and divide by n subjects
                %only includes first 45 rows to match matrix dim
                tempavg = eval(['epochs_eog_all.' conditions{ii}, '{1,1}(1:45,:)']);
            for iii = 2:length(sub_date.ID)
                tempavg = tempavg + eval(['epochs_eog_all.' conditions{ii}, '{iii,stim_index}(1:45,:)']);
            end
                tempavg = tempavg/length(sub_date.ID);

                epochs_eog_avgrast.(conditions{ii}){1, stim_index} = tempavg
            
        %For stim
        end
    
    %For conditions
    end  
    
%save(['../mat_data/epochs_eog_avgrast.mat'], 'epochs_eog_avgrast');

%% Carrier 60 plots

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
title({'Average EOG amplitude', 'after pulse only (n = 22)'});

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
title({'Average EOG amplitude with different', 'inter-stimulus interval (ISI) (n = 22)'});

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

title({'EOG amplitude after pulse only (n = 22)', '60dBA carrier'});

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

title({'EOG amplitude after gap-pulse (n = 22)', '60dBA carrier'});

xlabel('Inter-stimulus interval (ms)');
ylabel('Amplitude (uV)');

hold off

%% Carrier 70 plots

%Plot of mean EOG002 amplitude for all subjects PO70-stim
%NB "times" (X axis) defined manually, from C_process_timelockeds.m
figure
hold on
plot(mean(cell2mat(epochs_eog.PO70(:,1))));
plot(mean(cell2mat(epochs_eog.PO70(:,2))));
plot(mean(cell2mat(epochs_eog.PO70(:,3))));
plot(mean(cell2mat(epochs_eog.PO70(:,4))));
plot(mean(cell2mat(epochs_eog.PO70(:,5))));

ax = gca;

plot([35 35], [ax.YLim(1), ax.YLim(2)], '--r');
plot([50 50], [ax.YLim(1), ax.YLim(2)], '--r');

legend('75', '80', '85', '90', '95');
title({'Average EOG amplitude', 'after pulse only (n = 22)'});

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

% Plot average GP70/ISI
%NB "times" (X axis) defined manually, from C_process_timelockeds.m
figure
hold on
plot(mean(cell2mat(epochs_eog.GP70(:,1))));
plot(mean(cell2mat(epochs_eog.GP70(:,2))));
plot(mean(cell2mat(epochs_eog.GP70(:,3))));
plot(mean(cell2mat(epochs_eog.GP70(:,4))));

%plot(mean(cell2mat(epochs_eog.GO(:,1))), '--r', 'LineWidth', 1);

ax = gca;

plot([35 35], [ax.YLim(1), ax.YLim(2)], '--r');
plot([50 50], [ax.YLim(1), ax.YLim(2)], '--r');

times = [-20:1:60]/200*1000;

ax = gca;
ax.XTickLabel = times(1:10:end);
ax.XLim = [0 80];

legend('0','60','120','240');
title({'Average EOG amplitude with different', 'inter-stimulus interval (ISI) (n = 22)'});

[hleg,att] = legend('show');
title(hleg, 'ISI (ms)');
hleg.Location = 'southwest';

hold off


%Boxplot of maximum response amplitude PO70
boxlabel = [repmat({'75'}, length(epochs_eog_resp.PO70(:,1)), 1),
    repmat({'80'}, length(epochs_eog_resp.PO70(:,2)), 1),
    repmat({'85'}, length(epochs_eog_resp.PO70(:,3)), 1),
    repmat({'90'}, length(epochs_eog_resp.PO70(:,4)), 1),
    repmat({'95'}, length(epochs_eog_resp.PO70(:,5)), 1)];

figure
hold on
boxplot([cell2mat(epochs_eog_resp.PO70(:,1)), 
    cell2mat(epochs_eog_resp.PO70(:,2)), 
    cell2mat(epochs_eog_resp.PO70(:,3)), 
    cell2mat(epochs_eog_resp.PO70(:,4)), 
    cell2mat(epochs_eog_resp.PO70(:,5))], boxlabel);

title({'EOG amplitude after pulse only (n = 22)', '70dBA carrier'});

xlabel('Pulse level (dBA)');
ylabel('Amplitude (uV)');

hold off


%Boxplot Boxplot of maximum response amplitude GP60
boxlabel = [repmat({'ISI 0'}, length(epochs_eog_resp.PO70(:,1)), 1),
    repmat({'ISI 60'}, length(epochs_eog_resp.PO70(:,2)), 1),
    repmat({'ISI 120'}, length(epochs_eog_resp.PO70(:,3)), 1),
    repmat({'ISI 240'}, length(epochs_eog_resp.PO70(:,4)), 1)];

figure
hold on
boxplot([cell2mat(epochs_eog_resp.PO70(:,1)), 
    cell2mat(epochs_eog_resp.PO70(:,2)), 
    cell2mat(epochs_eog_resp.PO70(:,3)), 
    cell2mat(epochs_eog_resp.PO70(:,4))], boxlabel);

title({'EOG amplitude after gap-pulse (n = 5)', '70dBA carrier'});

xlabel('Inter-stimulus interval (ms)');
ylabel('Amplitude (uV)');

hold off

%% Habituation plots

EOG_threshold = 0.5*10^-4;

%PO60 lineplots - all trials for subjects, with positive threshold
for i = 1:length(sub_date.ID)

a = epochs_eog_all.PO60{i,2}(mean(epochs_eog_all.PO60{i,2}, 2) < EOG_threshold,:);

figure
hold on
    for j = 1:size(a,1)
        plot(a(j,:));
    end 
hold off

end

%GP60 lineplots - all trials for subjects, with positive threshold
for i = 1:length(sub_date.ID)

a = epochs_eog_all.GP60{i,2}(mean(epochs_eog_all.GP60{i,2}, 2) < EOG_threshold,:);

figure
hold on
    for j = 1:size(a,1)
        plot(a(j,:));
    end 
hold off

end

%PO60 raster for n = 22, with positive threshold
ha = tight_subplot(11,2,[.01 .05],[.1 .01],[.05 .01])
times = [-20:1:60]/200*1000;
          for ii = 1:22; axes(ha(ii));
              imagesc(epochs_eog_all.PO60{ii,6}(mean(epochs_eog_all.PO60{ii,6}, 2) < EOG_threshold,:))
              colormap bone;
              hold on
              plot([21 21], [0 50], 'k --');
              
              %should use "times" here, show x-axis if won't work?
              if ii == 21 | ii == 22
                ax = gca;
                ax.XTick = [1:10:80];
                ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];
              else
                ax = gca;
                ax.XTick = [1:10:80];
                ax.XTickLabel = [];
              end
              
              hold off
          end

%GP60 raster for n = 22, with positive threshold
ha = tight_subplot(11,2,[.01 .05],[.1 .01],[.05 .01])
times = [-20:1:60]/200*1000;
for ii = 1:22; axes(ha(ii));
  imagesc(epochs_eog_all.GP60{ii,3}(mean(epochs_eog_all.PO60{ii,3}, 2) < EOG_threshold,:))
  colormap bone;
  hold on
  plot([21 21], [0 50], 'k --');

  %should use "times" here
  ax = gca;
  ax.XTick = [1:10:80];
  ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];

  hold off
end
     
%The big ones
%times = [-20:1:60]/200*1000;

%PO60
figure('Position', [450 500 800 1000]);
ha = tight_subplot(6,1,[.01 .05],[.1 .01],[.05 .01])
for ii = 1:6; axes(ha(ii));
  imagesc(epochs_eog_avgrast.PO60{1,ii})
  colormap bone;
  hold on
  plot([21 21], [0 45], 'k --');

  %should use "times" here, show x-axis if won't work?
  if ii == 6
    ax = gca;
    ax.XTick = [1:10:80];
    ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];
  else
    ax = gca;
    ax.XTick = [1:10:80];
    ax.XTickLabel = [];
  end

  hold off
end


%GP60 + GO
figure('Position', [450 500 800 1000]);
ha = tight_subplot(6,1,[.01 .05],[.1 .01],[.05 .01])
for ii = 1:6; 
  if ii < 5
      axes(ha(ii));
      imagesc(epochs_eog_avgrast.GP60{1,ii})
      colormap bone;
      hold on
      plot([21 21], [0 45], 'k --');

      %should use "times" here, show x-axis if won't work?
      if ii == 4
        ax = gca;
        ax.XTick = [1:10:80];
        ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];
      else
        ax = gca;
        ax.XTick = [1:10:80];
        ax.XTickLabel = [];
      end
      
  elseif ii == 5
      axes(ha(ii));
      imagesc([])
      colormap bone;
      hold on

      %should use "times" here
      ax = gca;
      ax.XTick = [1:10:80];
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      
    elseif ii == 6
      axes(ha(ii));
      imagesc(epochs_eog_avgrast.GO{1,1})
      colormap bone;
      hold on
      plot([21 21], [0 45], 'k --');

      %should use "times" here
      ax = gca;
      ax.XTick = [1:10:80];
      ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];
 
  end
      
end

%PO70
figure('Position', [450 500 800 1000]);
ha = tight_subplot(5,1,[.025 .01],[.05 .05],[.05 .05])
for ii = 1:5; axes(ha(ii));
  imagesc(epochs_eog_avgrast.PO70{1,ii})
  colormap bone;
  hold on
  plot([21 21], [0 45], 'k --');

  %should use "times" here
  if ii == 5
    ax = gca;
    ax.XTick = [1:10:80];
    ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];
  else
    ax = gca;
    ax.XTick = [1:10:80];
    ax.XTickLabel = [];
  end

  hold off
end


%GP70 + GO
figure('Position', [450 500 800 1000]);
ha = tight_subplot(5,1,[.025 .01],[.05 .05],[.05 .05])
for ii = 1:5; 
  if ii < 5
      axes(ha(ii));
      imagesc(epochs_eog_avgrast.GP60{1,ii})
      colormap bone;
      hold on
      plot([21 21], [0 45], 'k --');

      %should use "times" here, show x-axis if won't work?
      if ii == 5
        ax = gca;
        ax.XTick = [1:10:80];
        ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];
      else
        ax = gca;
        ax.XTick = [1:10:80];
        ax.XTickLabel = [];
      end
      
    elseif ii == 5
      axes(ha(ii));
      imagesc(epochs_eog_avgrast.GO{1,2})
      colormap bone;
      hold on
      plot([21 21], [0 45], 'k --');

      %should use "times" here
      ax = gca;
      ax.XTick = [1:10:80];
      ax.XTickLabel = [-100 -50 0 50 100 150 200 250 300];
 
  end
      
end

