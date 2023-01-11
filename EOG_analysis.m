
load('../mat_data/timelockeds/epochs_eog_all.mat');
load(['../mat_data/timelockeds/epochs_eog_avgrast.mat']);
load('../mat_data/timelockeds/epochs_eog.mat');
load('../mat_data/timelockeds/epochs_eog_resp.mat');

% These two for Form9 manuscript
load('../mat_data/timelockeds/epochs_eog_all_clean_avg.mat');
load('../mat_data/timelockeds/epochs_eog_clean_resp.mat');

% Load one subject struct for time-info etc.
eog_timelockeds = load(['../mat_data/timelockeds/ID' num2str(sub_date.ID{1}) '/EOG/PO60_90_eog.mat'])
eog_timelockeds = eog_timelockeds.eog_timelockeds;

timevec = eog_timelockeds.time;

%% load EOG data to structure

epochs_eog = struct;

for i = 1:length(sub_date.ID)
    
    subinpath = ['../mat_data/timelockeds/ID' sub_date.ID{i} '/EOG/'];
    
    epochs_eog.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' conditions{ii} 'trig']));
    label = eval(['cond.' conditions{ii} 'label']);
    
        for stim_index = 1:nstim

        VAR = 'eog_timelockeds';
        T = load([subinpath char(label(stim_index)) '_eog' '.mat'], VAR);
        T = T.(VAR).avg;

        epochs_eog.(conditions{ii}){i, stim_index} = T
        
        %clear('T', 'VAR')
        
        %For stim
        end
    
    %For conditions
    end

%For subjects
end

%save(['../mat_data/timelockeds/epochs_eog.mat'], 'epochs_eog');

%% load EOG keeptrials data to structure

epochs_eog_all = struct;

for i = 1:length(sub_date.ID)
    
    subinpath = ['../mat_data/timelockeds/ID' sub_date.ID{i} '/EOG/'];
    
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
    
save(['../mat_data/timelockeds/epochs_eog_all.mat'], 'epochs_eog_all');

%% Clean EOG data
%Adapted from: https://se.mathworks.com/matlabcentral/answers/385226-how-to-use-the-data-brush-tool-to-automatically-save-selected-points-in-multiple-line-plots

% Pack up data
t = (1:165)'; %n of timepoints in data
names = regexp(cellstr(sprintf('trial_n_%d ',1:50)),' ','split');
names = names{1,1};

epochs_eog_all_clean = struct();

for i = 1:numel(conditions)

nstim = length(eval(['cond.' conditions{i} 'trig']));

    for ii = 1:nstim

        for iii = 1:numel(sub_date.ID)

            %load eog-trials in temp variable "y"
            y = epochs_eog_all.(conditions{i}){iii,ii}';
            
            % Plot data
            figure('Position',[800 800 1000 800])
            hLines = plot(t,y, 'k');
            xline([101 101]);
            xlim([1 165]);
            xticks([1:20:165]);
            xticklabels([-500:100:320]);
            
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
                    data = [t(ptsSelected) y(ptsSelected,k)];
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
            
            %Save plot of data with removed trials
            figure('Position',[800 800 1000 1000]); hold on;
            plot(t,y, 'k');
            if exist('colnum','var') == 1 %if any removed, plot them as red
                plot(t,y(:,[colnum]), 'r');
            end;
            xline([101 101]);
            xlim([1 165]);
            xticks([1:20:165]);
            xticklabels([-500:100:320]);

            saveas(gcf, ['../Analysis Output/EOG_cleaned/' sub_date.ID{iii} '_' conditions{i} '_' num2str(ii) '.jpg']);
            close;
            
            %if any brushed trials, remove from temp variable y
            if exist('colnum','var') == 1 %if any removed, plot them as red
                y(:,[colnum]) = [];
            end
    
            % Write y to new struct
            epochs_eog_all_clean.(conditions{i}){iii,ii} = y;

            clear('y');
        
        %for subjects
        end
    
    %for stim
    end

%for conditions
end

%save('../mat_data/timelockeds/epochs_eog_all_clean.mat', "epochs_eog_all_clean");




%% Plot cleaned EOG data

%load('../mat_data/timelockeds/epochs_eog_all_clean.mat');

for i = 1:numel(sub_date.ID)

    for ii = 1:numel(conditions)

        figure('Position', [500 500 1200 1100]); hold on;
        subplot(2,2,1)
        for iii = 1:numel(cond.([conditions{ii} 'label']))
            plot(epochs_eog_all.(conditions{ii}){i,iii}')
        end
        xline([101 101]);
        xticks([1:20:165]);
        xticklabels([-500:100:320]);
        xlim([1 165]);
        ylim([-6*10^-4 4*10^-4])
        title(['All trials ' (conditions{ii})]);
        
        subplot(2,2,2)
        for iii = 1:numel(cond.([conditions{ii} 'label']))
            plot(mean(epochs_eog_all.(conditions{ii}){i,iii}', 2))
        end
        xline([101 101]);
        xticks([1:20:165]);
        xticklabels([-500:100:320]);
        xlim([1 165]);
        ylim([-6*10^-4 4*10^-4])
        title(['Average per pulse level ' (conditions{ii})]);
        
        subplot(2,2,3)
        for iii = 1:numel(cond.([conditions{ii} 'label']))
            plot(epochs_eog_all_clean.(conditions{ii}){i,iii})
        end
        xline([101 101]);
        xticks([1:20:165]);
        xticklabels([-500:100:320]);
        xlim([1 165]);
        ylim([-6*10^-4 4*10^-4])
        title(['All trials ' (conditions{ii}) ' - cleaned']);
        
        subplot(2,2,4)
        for iii = 1:numel(cond.([conditions{ii} 'label']))
            plot(mean(epochs_eog_all_clean.(conditions{ii}){i,iii}, 2))
        end
        xline([101 101]);
        xticks([1:20:165]);
        xticklabels([-500:100:320]);
        xlim([1 165]);
        ylim([-6*10^-4 4*10^-4])
        title(['Average per pulse level ' (conditions{ii}) ' - cleaned']);

        %saveas(gcf, ['../Analysis Output/EOG_cleaned/' sub_date.ID{i} '_' conditions{ii} '.jpg']);
        close

    %for conditions
    end

%for subjects
end

%% Arrange cleaned data

%load('../mat_data/timelockeds/epochs_eog_all_clean.mat')

%Struct  to count how many trials were excluded
trials_left = struct();

%Struct for average of cleaned trials
epochs_eog_clean = struct();

for i = 1:numel(sub_date.ID)

    for ii = 1:numel(conditions)

        for iii = 1:numel(cond.([conditions{ii} 'label']))

            %Count n of excluded/kept trials
            trials_left.(conditions{ii}){i,iii} = size(epochs_eog_all_clean.(conditions{ii}){i,iii},2);

            %Collect mean for cleaned trials
            epochs_eog_clean.(conditions{ii}){i,iii} = mean(epochs_eog_all_clean.(conditions{ii}){i,iii}, 2)';

        end

    %for conditions
    end

%for subjects
end

%save('../mat_data/timelockeds/epochs_eog_all_clean_avg.mat', "epochs_eog_clean");

%% Find biggest response from structure

epochs_eog_resp = struct;

%Calculate max response for subjects
%NB, min between 50-150ms (sample 111-131)
for i = 1:length(sub_date.ID)

    epochs_eog_resp.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = numel(cond.([conditions{ii} 'trig']));
    label = cond.([conditions{ii} 'label']);
    
        for stim_index = 1:nstim
        
        %minresponse sample 111-131
        minresp = min(epochs_eog.(conditions{ii}){i,stim_index}(111:131));
        
        epochs_eog_resp.(conditions{ii})(i, stim_index) = minresp;
            
        %For stim
        end
    
    %For conditions
    end

%For subjects
end

%save(['../mat_data/timelockeds/epochs_eog_resp.mat'], 'epochs_eog_resp');

%Export EOG amps as CSV
%EOG_resp_csv = struct2table(epochs_eog_resp);
%writetable(EOG_resp_csv);

%% Find biggest response in cleaned data

epochs_eog_clean_resp = struct;

%From visual inspection of plots
toion = 111; %50ms
toioff = 147; %230ms

%Calculate max response for subjects
%NB, min between 50-150ms (sample 111-131)
for i = 1:numel(sub_date.ID)

    epochs_eog_clean_resp.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:numel(conditions)
    
    nstim = numel(cond.([conditions{ii} 'trig']));
    label = cond.([conditions{ii} 'label']);
    
        for stim_index = 1:nstim
        
        %minresponse in sample 111-131
        minresp = min(epochs_eog_clean.(conditions{ii}){i,stim_index}(toion:toioff));
        
        %replace "minresp" with [M, I] and use I for latency in separate
        %structure

        epochs_eog_clean_resp.(conditions{ii})(i, stim_index) = minresp;
            
        %For stim
        end
    
    %For conditions
    end

%For subjects
end

%save(['../mat_data/timelockeds/epochs_eog_clean_resp.mat'], 'epochs_eog_clean_resp');

%Export EOG amps as CSV
csvwrite('../R data/EOG_PO60', epochs_eog_clean_resp.PO60);
csvwrite('../R data/EOG_PO70', epochs_eog_clean_resp.PO70);
csvwrite('../R data/EOG_GP60', epochs_eog_clean_resp.GP60);
csvwrite('../R data/EOG_GP70', epochs_eog_clean_resp.GP70);

%% Find biggest response in all trials (conditions of interest)

%New empty structure
epochs_eog_allresp = struct;


%NB, min between 75-150ms (sample x-x)
%For all subjects
for i = 1:length(sub_date.ID)
    
    %List subject ID
    epochs_eog_allresp.subjects{i,1} = sub_date.ID{i};
    
    %For conditions 1 & 3: PO60 and GP60
    for ii = [1,3] 
        
        %Include only specific stim of interest (90dB PO and isi 240ms GP)
        if ii == 1; stim_index =  5; 
            elseif ii == 3; stim_index = 4;
        end;
        
        %Matrix of all subjects trials with stim
        minwin = cell2mat(epochs_eog_all.(conditions{ii})(i,stim_index));
        
        %Include first 45 trials for consitency
        ntrials = 45;
        
        %For each trial (row) take min between samples representing 75-150ms
        for iii = 1:ntrials

            %minresponse
            minresp = min(minwin(iii,35:50));
            epochs_eog_allresp.(conditions{ii}){i, iii} = minresp
        end
        
    %For conditions
    end

%For subjects
end

% Save as CSV
%csvwrite('../Analysis Output/EOG_response_PO60_90.csv', epochs_eog_allresp.PO60);
%csvwrite('../Analysis Output/EOG_response_GP60_i240.csv', epochs_eog_allresp.GP60);

%% Create struct of average rasters

epochs_eog_avgrast = struct;

for i = 1:length(sub_date.ID)

subinpath = ['../mat_data/timelockeds/ID' sub_date.ID{i} '/EOG/'];
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
    
end

%save(['../mat_data/timelockeds/epochs_eog_avgrast.mat'], 'epochs_eog_avgrast');

%% Average rasterplots

xtick = [1:10:165];
xticklab = [-500:50:320];

ytick = [0:10:45];
yticklab = [0:10:45];

triglinex = [101 101];
xrange = [13 161];

amp_ylim_low = -8*10^-5;
amp_ylim_high = 2*10^-5;

%Pulse Only average in rows
figure('Position', [100 500 2000 800]);
ha = tight_subplot(4,6,[.05 .015],[.2 .1],[.05 .05]);
colbar = [-8*10^-5 8*10^-5]; %set limits of color-gradient/"z-axis"

%PO70 raster row
for ii = 1:6; axes(ha(ii));
  
  if ii == 1
      axes(ha(ii));
      imagesc([])
      colormap parula;
      hold on
      
      xlim(xrange);
      
      ax = gca;
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      ax.XGrid = 'Off';
  
  elseif ii == 2
      imagesc(epochs_eog_avgrast.PO70{1,ii-1}, colbar)
      colormap parula;
      hold on
      plot(triglinex, [0 45], 'k --');
      
      xlim(xrange);
  
      ax = gca;
      
      ax.YTick = [ytick];
      ax.YTickLabel = [yticklab];         
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.XGrid = 'On';
      ax.FontSize = 12;
      
  elseif ii > 2
  
      imagesc(epochs_eog_avgrast.PO70{1,ii-1}, colbar)
      colormap parula;
      hold on
      plot(triglinex, [0 45], 'k --');

      xlim(xrange);

      ax = gca;
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      ax.XGrid = 'On';
      
  end
end

%PO70 amplitude line row
for ii = 1:6; axes(ha(ii+6));
      if ii == 1
      axes(ha(ii));
      imagesc([])
      colormap parula;
      hold on
      
      ax = gca;
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      ax.XGrid = 'On';
      
      elseif ii == 2
      plot(mean(epochs_eog_avgrast.PO70{1,ii-1}), 'k', 'LineWidth', 0.5);
      hold on
      plot(triglinex, [amp_ylim_low amp_ylim_high], 'k --');
      ylim([amp_ylim_low amp_ylim_high]);
      xlim(xrange);
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      ax.FontSize = 12;
      
      elseif ii > 2
      plot(mean(epochs_eog_avgrast.PO70{1,ii-1}), 'k', 'LineWidth', 0.5);
      hold on
      plot(triglinex, [amp_ylim_low amp_ylim_high], 'k --');
      ylim([amp_ylim_low amp_ylim_high]);
      xlim(xrange);
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      
      end
end

%PO60 raster in row
for ii = 1:6; axes(ha(ii+12));
  imagesc(epochs_eog_avgrast.PO60{1,ii}, colbar);
  colormap parula;
  hold on
  plot(triglinex, [0 45], 'k --');

  xlim(xrange);

  ax = gca;

  ax.YTick = [ytick];
  ax.YTickLabel = [yticklab];         
  ax.XTick = xtick;
  ax.XTickLabel = [];
  ax.XGrid = 'On';
  ax.FontSize = 12;
  
  if ii > 1
  plot(triglinex, [0 45], 'k --');

  xlim(xrange);

  ax = gca;

  ax.YTick = [ytick];
  ax.YTickLabel = [];         
  ax.XTick = xtick;
  ax.XTickLabel = [];
  ax.XGrid = 'On';
  
  end
  
end

%PO60 amplitude line row
for ii = 1:6; axes(ha(ii+18));
    
      plot(mean(epochs_eog_avgrast.PO60{1,ii}), 'k', 'LineWidth', 0.5);
      hold on
      
      plot(triglinex, [amp_ylim_low amp_ylim_high], 'k --');
      ylim([amp_ylim_low amp_ylim_high]);
      xlim(xrange);
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      
      if ii > 1
          ax = gca; 
          ax.YTickLabel = [];
      end
      
      ax = gca;
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      
      ax.XTick = [xtick];
      ax.XTickLabel = [xticklab];
      ax.XTickLabelRotation = 90;
      ax.FontSize = 12;

     ax.XTickLabel(2,:) = nan;
     ax.XTickLabel(4,:) = nan;
     ax.XTickLabel(6,:) = nan;
     ax.XTickLabel(8,:) = nan;
     ax.XTickLabel(10,:) = nan;
     ax.XTickLabel(12,:) = nan;
     ax.XTickLabel(14,:) = nan;
     ax.XTickLabel(16,:) = nan;

end


%Gap Only, Gap-Pulse in rows
figure('Position', [100 500 2000 800]);
ha = tight_subplot(4,6,[.05 .015],[.2 .1],[.05 .05]);
colbar = [-8*10^-5 8*10^-5]; %set limits of color-gradient/"z-axis"
%GP70 raster row
for ii = 1:6; axes(ha(ii));
    
      if ii == 1
      imagesc([]);
      colormap parula;
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.YTickLabel = [];
    
      elseif ii == 2

      imagesc(epochs_eog_avgrast.GO{1,2}, colbar) %NB manual cell reference
      colormap parula;
      hold on
      plot(triglinex, [0 45], 'k --');
      
      xlim(xrange);
  
      ax = gca;
      
      ax.YTick = [ytick];
      ax.YTickLabel = [yticklab];         
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.XGrid = 'On';
      ax.FontSize = 12;
      
      elseif ii > 2

      imagesc(epochs_eog_avgrast.GP70{1,ii-2}, colbar)
      colormap parula;
      hold on
      plot(triglinex, [0 45], 'k --');

      xlim(xrange);

      ax = gca;
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      ax.XGrid = 'On';
  end
      
end

%GP70 amplitude line row
for ii = 1:6; axes(ha(ii+6));
      if ii == 1
      axes(ha(ii+6));
      imagesc([])
      colormap parula;
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      
      hold on
      
      elseif ii == 2
          %static reference to GO cell
      plot(mean(epochs_eog_avgrast.GO{1,2}), 'k', 'LineWidth', 0.5);
      hold on
      plot(triglinex, [amp_ylim_low amp_ylim_high], 'k --');
      ylim([amp_ylim_low amp_ylim_high]);
      xlim(xrange);
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      ax.FontSize = 12;
      
      elseif ii > 2
      plot(mean(epochs_eog_avgrast.GP70{1,ii-2}), 'k', 'LineWidth', 0.5);
      hold on
      plot(triglinex, [amp_ylim_low amp_ylim_high], 'k --');
      ylim([amp_ylim_low amp_ylim_high]);
      xlim(xrange);
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      
      end
end

%GP60 raster row
for ii = 1:6; axes(ha(ii+12));
    
      if ii == 1
      axes(ha(ii+12));
      imagesc([])
      colormap parula;
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.YTickLabel = [];
    
      elseif ii == 2

      imagesc(epochs_eog_avgrast.GO{1,1}, colbar) %NB manual cell reference
      colormap parula;
      hold on
      plot(triglinex, [0 45], 'k --');
      
      xlim(xrange);
  
      ax = gca;
      
      ax.YTick = [ytick];
      ax.YTickLabel = [yticklab];         
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.XGrid = 'On';
      ax.FontSize = 12;
      
      elseif ii > 2

      imagesc(epochs_eog_avgrast.GP60{1,ii-2}, colbar)
      colormap parula;
      hold on
      plot(triglinex, [0 45], 'k --');
      
      xlim(xrange);
  
      ax = gca;
      
      ax.YTick = [ytick];
      ax.YTickLabel = [];         
      ax.XTick = xtick;
      ax.XTickLabel = [];
      ax.XGrid = 'On';
      ax.FontSize = 12;

  end
      
end

%GP60 amplitude line row
for ii = 1:6; axes(ha(ii+18));
      if ii == 1
      axes(ha(ii+18));
      imagesc([])
      colormap parula;
      
      ax = gca;
      ax.XTick = [xtick];
      ax.XTickLabel = [];
      ax.YTickLabel = [];
      
      hold on
      
      elseif ii == 2
          %static reference to GO cell
      plot(mean(epochs_eog_avgrast.GO{1,1}), 'k', 'LineWidth', 0.5);
      hold on
      plot(triglinex, [amp_ylim_low amp_ylim_high], 'k --');
      ylim([amp_ylim_low amp_ylim_high]);
      xlim(xrange);
      
      ax = gca;
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      
      ax.XTick = [xtick];
      ax.XTickLabel = [xticklab];
      ax.XTickLabelRotation = 90;
      ax.FontSize = 12;

      ax.XTickLabel(2,:) = nan;
      ax.XTickLabel(4,:) = nan;
      ax.XTickLabel(6,:) = nan;
      ax.XTickLabel(8,:) = nan;
      ax.XTickLabel(10,:) = nan;
      ax.XTickLabel(12,:) = nan;
      ax.XTickLabel(14,:) = nan;
      ax.XTickLabel(16,:) = nan;
      
      elseif ii > 2
      plot(mean(epochs_eog_avgrast.GP60{1,ii-2}), 'k', 'LineWidth', 0.5);
      hold on
      plot(triglinex, [amp_ylim_low amp_ylim_high], 'k --');
      ylim([amp_ylim_low amp_ylim_high]);
      xlim(xrange);
      
      ax = gca;
      ax.YTickLabel = [];
      ax.XGrid = 'On';
      ax.YGrid = 'On';
      
      ax.XTick = [xtick];
      ax.XTickLabel = [xticklab];
      ax.XTickLabelRotation = 90;
      ax.FontSize = 12;

      ax.XTickLabel(2,:) = nan;
      ax.XTickLabel(4,:) = nan;
      ax.XTickLabel(6,:) = nan;
      ax.XTickLabel(8,:) = nan;
      ax.XTickLabel(10,:) = nan;
      ax.XTickLabel(12,:) = nan;
      ax.XTickLabel(14,:) = nan;
      ax.XTickLabel(16,:) = nan;
      
      end
end

% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters PO_avg
%The first two lines measure the size of your figure (in inches). The next line configures the print paper size to fit the figure size. 
%The last line uses the print command and exports a vector pdf document as the output.




%% Permutetest of EOG

%create matrix (that works for permutest) of subject responses
for i = 1:numel(sub_date.ID)
    PO60avg(:,i) = epochs_eog.PO60{i,5}'; %NB condition 5 = 90dB pulse
    GP60avg(:,i) = epochs_eog.GP60{i,4}';
end

[clusters60, p_values60, t_sums60, permutation_distribution60] = permutest(PO60avg, GP60avg, true, 0.01, 10000, true, inf);

%create matrix (that works for permutest) of subject responses
for i = 1:numel(sub_date.ID)
    PO70avg(:,i) = epochs_eog.PO70{i,4}'; %NB condition 4 = 90dB pulse
    GP70avg(:,i) = epochs_eog.GP70{i,4}';
end

[clusters70, p_values70, t_sums70, permutation_distribution70] = permutest(PO70avg, GP70avg, true, 0.01, 10000, true, inf);

%Clean responses:
%create matrix (that works for permutest) of subject responses
for i = 1:numel(sub_date.ID)
    PO60_clean_avg(:,i) = epochs_eog_clean.PO60{i,5}'; %NB condition 5 = 90dB pulse
    GP60_clean_avg(:,i) = epochs_eog_clean.GP60{i,4}';
end

[clusters60c, p_values60c, t_sums60c, permutation_distribution60c] = permutest(PO60_clean_avg, GP60_clean_avg, true, 0.01, 10000, true, inf);

%create matrix (that works for permutest) of subject responses
for i = 1:numel(sub_date.ID)
    PO70_clean_avg(:,i) = epochs_eog.PO70{i,4}'; %NB condition 4 = 90dB pulse
    GP70_clean_avg(:,i) = epochs_eog.GP70{i,4}';
end

[clusters70c, p_values70c, t_sums70c, permutation_distribution70c] = permutest(PO70_clean_avg, GP70_clean_avg, true, 0.01, 10000, true, inf);

%% Manuscript plots_v1

%Make a grand average EOG response
for i = 1:numel(conditions)
    for ii = 1:numel(cond.([(conditions{i}) 'label']))
        for iii = 1:numel(sub_date.ID)
            temp(iii,:) = epochs_eog.(conditions{i}){iii,ii};
        end
    %SEM calculation goes here
    EOG_GA.(conditions{i})(ii,:) = mean(temp,1);
    clear temp
    end
end

% Inherited from raster plots
xtick = [1:20:165];
xticklab = [-500:100:320];
triglinex = [101 101];
xrange = [41 161];

txtsize = 10;

lineylims = [-8*10^-5 1*10^-5];
boxylims = [-2*10^-4 0.5*10^-4];

%Tableau medium 10 palette
colors = [173 139 201;
    168 120 110;
    114 158 206;
    255 158 74;
    237 102 93;
    103 191 92;
    237 151 202;
    162 162 162;
    205 204 93;
    109 204 218]/256;

%No idea why the colors are red backwards for boxplot
boxcolors = flip(colors(1:6,:));
isiboxcolors = flip(colors(7:10,:));

figure('Units', 'centimeters', 'Position',  [5 5 25 30]); subplot(4,3,1); 
hold on;

%PO60
for i = 1:6
    plot(EOG_GA.PO60(i,:), 'Color', colors(i,:))
end

plot(triglinex, lineylims, 'k --');

ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [];
ax.XGrid = 'on';
ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);

ylabel({"60dB carrier", "EOG Amplitude (V)"})

legend({"70", "75", "80", "85", "90", "95"}, 'Location', 'southwest');
legend('boxoff');

subplot(4,3,2); boxplot(epochs_eog_resp.PO60(:, :), [70:5:95], 'Symbol', 'ok');
ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
xlabel("Pulse level")

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:6
    patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

%Inhib PO6090 & ISI240
subplot(4,3,3); hold on;
plot(mean(PO60avg,2), 'Color', colors(5,:)); %Manually specify colors to match other plots
plot(mean(GP60avg,2), 'Color', colors(10,:)); % -''-
plot(clusters60{:}, repmat(0.5*10^-5, 1, length(clusters60{:})), 'red', 'LineWidth', 3)

ylim(lineylims);
xlim(xrange);

h = get(gca,'Children');
set(gca,'Children',[h(3) h(2) h(1)])

ax = gca;
ax.FontSize = txtsize;
ax.XTick = xtick;
ax.XTickLabel = xticklab;
ax.XGrid = 'on';
ax.XTickLabelRotation = 90;

xlabel("Time (ms)");

legend({"Cluster, p = 0.01", "Gap + Pulse (ISI 240 ms)", "Pulse (90 dB)"}, 'Location', 'southwest');
legend('boxoff');


%PO70
subplot(4,3,4); hold on;
for i = 1:5
    plot(EOG_GA.PO70(i,:), 'Color', colors(i+1,:)) %+1 to color to match PO60
end
plot(triglinex, lineylims, 'k --');
ylim(lineylims)

ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [];
ax.XGrid = 'on';
ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);
ylabel({"70dB carrier", "EOG Amplitude (V)"})

%NaNs to pad missing 70dB pulse in 70dB carrier
subplot(4,3,5); boxplot([repmat(NaN, 1, 22)' epochs_eog_resp.PO70(:, 1:5)], [70:5:95], 'Symbol', 'ok');
ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
xlabel("Pulse level")

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:5
    patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

%GP60
subplot(4,3,7); hold on;
for i = 1:4
    plot(EOG_GA.GP60(i,:), 'Color', colors(i+6,:))
end
plot(triglinex, lineylims, 'k --');
ylim(lineylims)

ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [];
ax.XGrid = 'on';
ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);
ylabel({"60dB carrier", "EOG Amplitude (V)"})

legend({"0 ms", "60 ms", "120 ms", "240 ms"}, 'Location', 'southwest');
legend('boxoff');

subplot(4,3,8); boxplot(epochs_eog_resp.GP60(:, :), [0 60 120 240], 'Symbol', 'ok');
ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
xlabel("ISI duration")

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:4
    patch(get(h(j),'XData'),get(h(j),'YData'),isiboxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

%GP70
subplot(4,3,10); hold on;
for i = 1:4
    plot(EOG_GA.GP70(i,:), 'Color', colors(i+6,:))
end
plot(triglinex, lineylims, 'k --');
ylim(lineylims)

ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [xticklab];
ax.XGrid = 'on';
ax.XTickLabelRotation = 90;

ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);
ylabel({"70dB carrier", "EOG Amplitude (V)"})
xlabel("Time (ms)")

subplot(4,3,11); boxplot(epochs_eog_resp.GP70(:, :), [0 60 120 240], 'Symbol', 'ok');
ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
xlabel("ISI duration")

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:4
    patch(get(h(j),'XData'),get(h(j),'YData'),isiboxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

%saveas(gcf, ['../Analysis Output/MEG_manus_EOG.svg']);

%% Mauscript plots_V2 60 dB carrier

% Form 9, fig 1
% V2: Using cleaned up EOG responses

%load('../mat_data/timelockeds/epochs_eog_all_clean_avg.mat');

%Make a grand average EOG response
for i = 1:numel(conditions)
    for ii = 1:numel(cond.([(conditions{i}) 'label']))
        for iii = 1:numel(sub_date.ID)
            temp(iii,:) = epochs_eog_clean.(conditions{i}){iii,ii};
        end
    %SEM calculation goes here
    EOG_GA_clean.(conditions{i})(ii,:) = mean(temp,1);
    clear temp
    end
end

% Inherited from raster plots
xtick = [1:20:165];
xticklab = [-500:100:320];
triglinex = [101 101];
xrange = [41 161];

txtsize = 12;
linew = 1.5;

lineylims = [-9*10^-5 1*10^-5];
boxylims = [-2*10^-4 0.5*10^-4];

toion = 111; %50ms
toioff = 147; %230ms

%Tableau medium 10 palette
colors = [173 139 201;
    168 120 110;
    114 158 206;
    255 158 74;
    237 102 93;
    103 191 92;
    237 151 202;
    162 162 162;
    205 204 93;
    109 204 218]/256;

%No idea why the colors are read backwards for boxplot
boxcolors = flip(colors(1:6,:));
isiboxcolors = flip(colors(7:10,:));

figure('Units', 'centimeters', 'Position',  [5 5 30 20], 'Renderer','painters');
tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

%PO60
nexttile; hold on;
patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:6
    plot(EOG_GA_clean.PO60(i,:), 'Color', colors(i,:), 'LineWidth', linew)
end

plot(triglinex, lineylims, 'k --');
ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [xticklab];
ax.XGrid = 'on';
ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);
xlabel("Time (ms)");

ylabel({"60dB carrier", "EOG Amplitude (V)"})

title({'Pulse only trials', 'Different pulse level'});

legend({"", "70", "75", "80", "85", "90", "95"}, 'Location', 'southwest'); %First one empty to skip patch
legend('boxoff');

nexttile; hold on;
boxplot(epochs_eog_clean_resp.PO60(:, :), [70:5:95], 'Symbol', 'ok');

title({'Pulse only response amplitude', 'Different pulse level'});

ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
ax.Box = 'off';

set(findobj(gca,'type','line'),'lineStyle','-');

xlabel("Pulse level");

%[normx, normy] = coord2norm(ax, [5 5], [-18*10^-5 -20*10^-5]); %https://se.mathworks.com/matlabcentral/fileexchange/54254-coord2norm
%annotation('arrow', normx, normy);

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:6
    patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

title({'Pulse only response amplitude', 'Different pulse level'});



%GP60
nexttile; hold on;
patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:4
    plot(EOG_GA_clean.GP60(i,:), 'Color', colors(i+6,:), 'LineWidth', linew)
end
plot(triglinex, lineylims, 'k --');
ylim(lineylims)

ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [xticklab];
ax.XGrid = 'on';
ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);

xlabel("Time (ms)");
ylabel({"60dB carrier", "EOG Amplitude (V)"})

legend({"", "0 ms", "60 ms", "120 ms", "240 ms"}, 'Location', 'southwest'); %First one empty to skip patch
legend('boxoff');

title({'Gap + Pulse trials', 'Different ISI'});

nexttile; boxplot(epochs_eog_clean_resp.GP60(:, :), [0 60 120 240], 'Symbol', 'ok');
ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
xlabel("ISI duration")
ax.Box = 'off';

set(findobj(gca,'type','line'),'lineStyle','-');

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:4
    patch(get(h(j),'XData'),get(h(j),'YData'),isiboxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

title({'Gap + Pulse response amplitude', 'Different ISI'});


% %PO70
% nexttile; hold on;
% patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
% for i = 1:5
%     plot(EOG_GA_clean.PO70(i,:), 'Color', colors(i+1,:), 'LineWidth', linew) %+1 to color to match PO60
% end
% plot(triglinex, lineylims, 'k --');
% ylim(lineylims)
% 
% ax = gca;
% ax.XTick = xtick;
% ax.XTickLabel = [xticklab];
% ax.XGrid = 'on';
% ax.FontSize = txtsize;
% ylim(lineylims);
% xlim(xrange);
% ylabel({"70dB carrier", "EOG Amplitude (V)"})
% xlabel("Time (ms)");
% 
% legend({"", "75", "80", "85", "90", "95"}, 'Location', 'southwest'); %First one empty to skip patch
% legend('boxoff');
% 
% %NaNs to pad missing 70dB pulse in 70dB carrier
% nexttile; boxplot([repmat(NaN, 1, 22)' epochs_eog_clean_resp.PO70(:, 1:5)], [70:5:95], 'Symbol', 'ok');
% ylim(boxylims)
% ax = gca;
% ax.FontSize = txtsize;
% xlabel("Pulse level")
% ax.Box = 'off';
% 
% set(findobj(gca,'type','line'),'lineStyle','-');
% 
% h = findobj(gcf,'tag','Outliers');
% set(h,'MarkerSize',4);
% 
% h = findobj(gca,'Tag','Box');
% for j=1:5
%     patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',.5);
% end
% 
% lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% set(lines, 'Color', 'k');
% 
% %GP70
% nexttile; hold on;
% patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
% for i = 1:4
%     plot(EOG_GA_clean.GP70(i,:), 'Color', colors(i+6,:), 'LineWidth', linew)
% end
% plot(triglinex, lineylims, 'k --');
% ylim(lineylims)
% 
% ax = gca;
% ax.XTick = xtick;
% ax.XTickLabel = [xticklab];
% ax.XGrid = 'on';
% 
% ax.FontSize = txtsize;
% ylim(lineylims);
% xlim(xrange);
% xlabel("Time (ms)");
% 
% legend({"", "0 ms", "60 ms", "120 ms", "240 ms"}, 'Location', 'southwest'); %First one empty to skip patch
% legend('boxoff');
% 
% nexttile; boxplot(epochs_eog_clean_resp.GP70(:, :), [0 60 120 240], 'Symbol', 'ok');
% ylim(boxylims)
% ax = gca;
% ax.FontSize = txtsize;
% xlabel("ISI duration")
% ax.Box = 'off';
% 
% set(findobj(gca,'type','line'),'lineStyle','-');
% 
% h = findobj(gcf,'tag','Outliers');
% set(h,'MarkerSize',4);
% 
% h = findobj(gca,'Tag','Box');
% for j=1:4
%     patch(get(h(j),'XData'),get(h(j),'YData'),isiboxcolors(j,:),'FaceAlpha',.5);
% end
% 
% lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% set(lines, 'Color', 'k');

%Inhib PO6090 & ISI240
% nexttile; hold on;
% 
% patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
% plot(EOG_GA_clean.PO60(5,:), 'Color', colors(5,:), 'LineWidth', linew); %Manually specify colors to match other plots
% plot(EOG_GA_clean.GP60(4,:), 'Color', colors(10,:), 'LineWidth', linew); % -''-
% 
% % %plot clusters
% % for j = 1:numel(clusters60c)
% %     plot(clusters60c{j}, repmat(0.5*10^-5, 1, numel(clusters60c{j})), 'red', 'LineWidth', 3)
% % end
% 
% ylim(lineylims);
% xlim(xrange);
% 
% ax = gca;
% ax.FontSize = txtsize;
% ax.XTick = xtick;
% ax.XTickLabel = xticklab;
% ax.XGrid = 'on';
% xline([101], '--k')
% 
% xlabel("Time (ms)");
% 
% legend({"", "Pulse (90 dB)", "Gap + Pulse (ISI 240 ms)"}, 'Location', 'southwest');
% legend('boxoff');
% 
% title({'Inhibition of EOG response', '60 dB carrier'});

%%%%%%%%%%%%%
% %Inhib PO7090 & ISI240
% nexttile; hold on;
% plot(mean(PO70_clean_avg,2), 'Color', colors(5,:)); %Manually specify colors to match other plots
% plot(mean(GP70_clean_avg,2), 'Color', colors(10,:)); % -''-
% 
% %plot clusters
% for j = 1:numel(clusters70c)
%     plot(clusters70c{j}, repmat(0.5*10^-5, 1, numel(clusters70c{j})), 'red', 'LineWidth', 3)
% end
% 
% ylim(lineylims);
% xlim(xrange);
% 
% ax = gca;
% ax.FontSize = txtsize;
% ax.XTick = xtick;
% ax.XTickLabel = xticklab;
% ax.XGrid = 'on';
% 
% xlabel("Time (ms)");
% 
% legend({"Pulse (90 dB)", "Gap + Pulse (ISI 240 ms)", "Cluster, p = 0.01"}, 'Location', 'southwest');
% legend('boxoff');
% 
% title({'Inhibition of EOG response', '70 dB carrier'});

saveas(gcf, ['../Analysis Output/Fig1_form9.svg']);

%% Mauscript plots_V2 70 dB carrier

% Form 9, fig 1
% V2: Using cleaned up EOG responses

%load('../mat_data/timelockeds/epochs_eog_all_clean_avg.mat');

%Make a grand average EOG response
for i = 1:numel(conditions)
    for ii = 1:numel(cond.([(conditions{i}) 'label']))
        for iii = 1:numel(sub_date.ID)
            temp(iii,:) = epochs_eog_clean.(conditions{i}){iii,ii};
        end
    %SEM calculation goes here
    EOG_GA_clean.(conditions{i})(ii,:) = mean(temp,1);
    clear temp
    end
end

% Inherited from raster plots
xtick = [1:20:165];
xticklab = [-500:100:320];
triglinex = [101 101];
xrange = [41 161];

txtsize = 12;
linew = 1.5;

lineylims = [-9*10^-5 1*10^-5];
boxylims = [-2*10^-4 0.5*10^-4];

toion = 111; %50ms
toioff = 147; %230ms

%Tableau medium 10 palette
colors = [173 139 201;
    168 120 110;
    114 158 206;
    255 158 74;
    237 102 93;
    103 191 92;
    237 151 202;
    162 162 162;
    205 204 93;
    109 204 218]/256;

%No idea why the colors are read backwards for boxplot
boxcolors = flip(colors(1:6,:));
isiboxcolors = flip(colors(7:10,:));

figure('Units', 'centimeters', 'Position',  [5 5 30 20], 'Renderer','painters');
tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');


%PO70
nexttile; hold on;
patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:5
    plot(EOG_GA_clean.PO70(i,:), 'Color', colors(i+1,:), 'LineWidth', linew) %+1 to color to match PO60
end
plot(triglinex, lineylims, 'k --');
ylim(lineylims)

ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [xticklab];
ax.XGrid = 'on';
ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);
ylabel({"70dB carrier", "EOG Amplitude (V)"})
xlabel("Time (ms)");

legend({"", "75", "80", "85", "90", "95"}, 'Location', 'southwest'); %First one empty to skip patch
legend('boxoff');

title({'Pulse only trials', 'Different pulse level'});

%NaNs to pad missing 70dB pulse in 70dB carrier
nexttile; boxplot(epochs_eog_clean_resp.PO70(:, 1:5), [75:5:95], 'Symbol', 'ok');
ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
xlabel("Pulse level")
ax.Box = 'off';

set(findobj(gca,'type','line'),'lineStyle','-');

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:5
    patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

title({'Pulse only response amplitude', 'Different pulse level'});

%GP70
nexttile; hold on;
patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:4
    plot(EOG_GA_clean.GP70(i,:), 'Color', colors(i+6,:), 'LineWidth', linew)
end
plot(triglinex, lineylims, 'k --');
ylim(lineylims)

ax = gca;
ax.XTick = xtick;
ax.XTickLabel = [xticklab];
ax.XGrid = 'on';

ax.FontSize = txtsize;
ylim(lineylims);
xlim(xrange);
xlabel("Time (ms)");
ylabel({"70dB carrier", "EOG Amplitude (V)"})

legend({"", "0 ms", "60 ms", "120 ms", "240 ms"}, 'Location', 'southwest'); %First one empty to skip patch
legend('boxoff');

title({'Gap + Pulse trials amplitude', 'Different ISI'});

nexttile; boxplot(epochs_eog_clean_resp.GP70(:, :), [0 60 120 240], 'Symbol', 'ok');
ylim(boxylims)
ax = gca;
ax.FontSize = txtsize;
xlabel("ISI duration")
ax.Box = 'off';

set(findobj(gca,'type','line'),'lineStyle','-');

h = findobj(gcf,'tag','Outliers');
set(h,'MarkerSize',4);

h = findobj(gca,'Tag','Box');
for j=1:4
    patch(get(h(j),'XData'),get(h(j),'YData'),isiboxcolors(j,:),'FaceAlpha',.5);
end

lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

title({'Gap + Pulse response amplitude', 'Different ISI'});

%Inhib PO6090 & ISI240
% nexttile; hold on;
% 
% patch('Faces', [1 2 3 4], 'Vertices', [toion lineylims(1); toion lineylims(2); toioff lineylims(2); toioff lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
% plot(EOG_GA_clean.PO60(5,:), 'Color', colors(5,:), 'LineWidth', linew); %Manually specify colors to match other plots
% plot(EOG_GA_clean.GP60(4,:), 'Color', colors(10,:), 'LineWidth', linew); % -''-
% 
% % %plot clusters
% % for j = 1:numel(clusters60c)
% %     plot(clusters60c{j}, repmat(0.5*10^-5, 1, numel(clusters60c{j})), 'red', 'LineWidth', 3)
% % end
% 
% ylim(lineylims);
% xlim(xrange);
% 
% ax = gca;
% ax.FontSize = txtsize;
% ax.XTick = xtick;
% ax.XTickLabel = xticklab;
% ax.XGrid = 'on';
% xline([101], '--k')
% 
% xlabel("Time (ms)");
% 
% legend({"", "Pulse (90 dB)", "Gap + Pulse (ISI 240 ms)"}, 'Location', 'southwest');
% legend('boxoff');
% 
% title({'Inhibition of EOG response', '60 dB carrier'});

%%%%%%%%%%%%%
% %Inhib PO7090 & ISI240
% nexttile; hold on;
% plot(mean(PO70_clean_avg,2), 'Color', colors(5,:)); %Manually specify colors to match other plots
% plot(mean(GP70_clean_avg,2), 'Color', colors(10,:)); % -''-
% 
% %plot clusters
% for j = 1:numel(clusters70c)
%     plot(clusters70c{j}, repmat(0.5*10^-5, 1, numel(clusters70c{j})), 'red', 'LineWidth', 3)
% end
% 
% ylim(lineylims);
% xlim(xrange);
% 
% ax = gca;
% ax.FontSize = txtsize;
% ax.XTick = xtick;
% ax.XTickLabel = xticklab;
% ax.XGrid = 'on';
% 
% xlabel("Time (ms)");
% 
% legend({"Pulse (90 dB)", "Gap + Pulse (ISI 240 ms)", "Cluster, p = 0.01"}, 'Location', 'southwest');
% legend('boxoff');
% 
% title({'Inhibition of EOG response', '70 dB carrier'});

saveas(gcf, ['../Analysis Output/FigSup1_form9.svg']);

%% Supplementary figure of latency variability for all subjects all stimuli

%load('../mat_data/timelockeds/epochs_eog_all_clean.mat');

timevec = -500:5:320; %Vector of timepoints in ms

xrange = [1 165];
ylims = [-6*10^-4 2*10^-4];

trigidx = find(timevec == 0); %Trigger index, i.e. time zero

xtick = [1:20:165];
xticklab = [-500:100:320];
triglinex = [trigidx trigidx];

for i = 1:numel(conditions)

    for ii = 1:numel(cond.([(conditions{i}) 'label']))
        
        figure('Position', [800 800 1000 1100], 'Renderer', 'painters');
        tiledlayout(6,4, 'TileSpacing','compact', 'Padding','compact'); 

        for iii = 1:numel(sub_date.ID)

            tempdat = epochs_eog_all_clean.(conditions{i}){iii,ii};
            
            SDvec = std(tempdat, 0, 2)';
            meanvec = mean(tempdat, 2)';
            
            SDposamp = meanvec + SDvec;
            SDnegamp = meanvec - SDvec;
            
            x = 1:numel(meanvec);
            x2 = [x, fliplr(x)];
            inBetween = [SDnegamp, fliplr(SDposamp)]';
            
            %find minimum magnitude and index between trigidx (time zero) and end
            [M,I] = min(tempdat(trigidx:numel(timevec),:),[], 1);
            
            %mean latency as ms
            meanlat = mean(timevec(I+trigidx-1));
            meanlatidx = interp1(timevec,1:numel(timevec),meanlat,'nearest');
            
            % Latency standard deviation
            SDlat = std(timevec(I+trigidx-1));
            
            % Index of mean +/- SD
            SDposlat = interp1(timevec,1:numel(timevec),(meanlat + SDlat),'nearest');
            SDneglat = interp1(timevec,1:numel(timevec),(meanlat - SDlat),'nearest');
            
            % Index of mean +/- SD
            SD_2_poslat = interp1(timevec,1:numel(timevec),(meanlat + SDlat*2),'nearest');
            SD_2_neglat = interp1(timevec,1:numel(timevec),(meanlat - SDlat*2),'nearest');           
            
            nexttile; hold on;

            %Latency patch + text of 2SD range
            %patch([SD_2_neglat SD_2_neglat SD_2_poslat SD_2_poslat], [ylims(1) ylims(2) ylims(2) ylims(1)], [0 0 1], 'FaceAlpha', 0.1, 'LineStyle', 'none');
            patch([SDneglat SDneglat SDposlat SDposlat], [ylims(1) ylims(2) ylims(2) ylims(1)], [0 0 1], 'FaceAlpha', 0.1, 'LineStyle', 'none');
            
            text(5, -5.5*10^-4, ['2 SD range: ' num2str((round(meanlat - SDlat*2,0))) ' - ' num2str((round(meanlat + SDlat*2,0)))], 'HorizontalAlignment', 'left', 'Color', 'blue', 'FontSize', 8);
            
            %Amplitude SD fills
            fill(x2, inBetween*2, [1 0 0], 'FaceAlpha', 0.15, 'LineStyle', 'none');
            fill(x2, inBetween, [1 0 0], 'FaceAlpha', 0.25, 'LineStyle', 'none');
            
            %Trials
            plot(tempdat, 'Color', [0 0 0 0.5])
            
            %Min amp dot
            plot(I+trigidx-1,M, 'r.')
            
            %Mean of trials
            plot(meanvec, 'r', 'LineWidth', 2)
            
            %Trigger marker
            plot(triglinex, ylims, 'k --');
            
            %Mean latency
            plot([meanlatidx meanlatidx], ylims, 'b --');
            
            xlim(xrange);
            ylim(ylims);
            
            ax = gca;
            ax.XTick = xtick;
            ax.XTickLabel = xticklab;

            title(sub_date.ID{iii})

            clear tempdat;

        end %subject
    
        saveas(gcf, ['../Analysis Output/EOG_variance/' cond.([(conditions{i}) 'label']){ii} '.svg'])
        close;
    end %stim

end %conditions

%% EOG Raster plot for Chris

if exist('epochs_eog_all_clean', 'var') == 0
    load('../mat_data/timelockeds/epochs_eog_all_clean.mat');
end

%index of max for PO60_90
if exist('epochs_eog_clean_resp', 'var') == 0
    load('../mat_data/timelockeds/epochs_eog_clean_resp.mat');
end

%Plot
colbar = [-1.5*10^-4 0.3*10^-4];
ylimamp = [-3*10^-4 1*10^-4];

figure('Position', [600 200 1400 1200]); tiledlayout(5,2, 'TileSpacing', 'compact');

%PO60_90
nexttile; hold on;
for i = 1:22
    temp(i, 1:165) = mean(epochs_eog_all_clean.PO60{i,5}, 2);
    plot(mean(epochs_eog_all_clean.PO60{i,5}, 2), 'Color', [0 0 0 0.5])
end
plot(mean(temp', 2), 'Color', [1 0 0], 'LineWidth', 2)

xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
xline(101)
ylim(ylimamp)
title('EOG 90dB Pulse only')
ylabel('EOG amp (uV)')
set(gca, 'XGrid', 'on');

clear ind
%Sort according to mean amplitud 75-150ms and keep index of sort
[val, ind] = sort(mean(temp(:,116:131), 2), 'ascend');
temp = temp(ind,:);

nexttile;
imagesc(temp, colbar); hold on;
xline(101)
colormap(viridis)
ylim([1 22])
xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
title('EOG 90dB Pulse only')
ylabel('subjects (n = 22)')
set(gca, 'XGrid', 'on');
colorbar;


%ISI_240
nexttile; hold on;
for i = 1:22
    temp(i, 1:165) = mean(epochs_eog_all_clean.GP60{i,4}, 2);
    plot(mean(epochs_eog_all_clean.GP60{i,4}, 2), 'Color', [0 0 0 0.5])
end
plot(mean(temp', 2), 'Color', [1 0 0], 'LineWidth', 2)

xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
xline(101)
ylim(ylimamp)
title('EOG ISI 240ms')
ylabel('EOG amp (uV)')
set(gca, 'XGrid', 'on');

temp = temp(ind,:);

nexttile;
imagesc(temp, colbar); hold on;
xline(101)
colormap(viridis)
ylim([1 22])
xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
title('EOG ISI 240ms')
ylabel('subjects (n = 22)')
set(gca, 'XGrid', 'on');

%ISI_120
nexttile; hold on;
for i = 1:22
    temp(i, 1:165) = mean(epochs_eog_all_clean.GP60{i,3}, 2);
    plot(mean(epochs_eog_all_clean.GP60{i,3}, 2), 'Color', [0 0 0 0.5])
end
plot(mean(temp', 2), 'Color', [1 0 0], 'LineWidth', 2)

xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
xline(101)
ylim(ylimamp)
title('EOG ISI 120ms')
ylabel('EOG amp (uV)')
set(gca, 'XGrid', 'on');

temp = temp(ind,:);

nexttile;
imagesc(temp, colbar); hold on;
xline(101)
colormap(viridis)
ylim([1 22])
xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
title('EOG ISI 120ms')
ylabel('subjects (n = 22)')
set(gca, 'XGrid', 'on');

%ISI_60
nexttile; hold on;
for i = 1:22
    temp(i, 1:165) = mean(epochs_eog_all_clean.GP60{i,2}, 2);
    plot(mean(epochs_eog_all_clean.GP60{i,2}, 2), 'Color', [0 0 0 0.5])
end
plot(mean(temp', 2), 'Color', [1 0 0], 'LineWidth', 2)

xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
xline(101)
ylim(ylimamp)
title('EOG ISI 60ms')
ylabel('EOG amp (uV)')
set(gca, 'XGrid', 'on');

temp = temp(ind,:);

nexttile;
imagesc(temp, colbar); hold on;
xline(101)
colormap(viridis)
ylim([1 22])
xlim([1 165])
xticks([1:20:165]);
xticklabels([]);
title('EOG ISI 60ms')
ylabel('subjects (n = 22)')
set(gca, 'XGrid', 'on');

%ISI_0
nexttile; hold on;
for i = 1:22
    temp(i, 1:165) = mean(epochs_eog_all_clean.GP60{i,1}, 2);
    plot(mean(epochs_eog_all_clean.GP60{i,1}, 2), 'Color', [0 0 0 0.5])
end
plot(mean(temp', 2), 'Color', [1 0 0], 'LineWidth', 2)

xlim([1 165])
xticks([1:20:165]);
xticklabels([-500:100:320]);
xline(101)
ylim(ylimamp)
title('EOG ISI 0ms')
ylabel('EOG amp (uV)')
xlabel('Time (ms)')
set(gca, 'XGrid', 'on');

temp = temp(ind,:);

nexttile;
imagesc(temp, colbar); hold on;
xline(101)
colormap(viridis)
ylim([1 22])
xlim([1 165])
xticks([1:20:165]);
xticklabels([-500:100:320]);
title('EOG ISI 0ms')
ylabel('subjects (n = 22)')
xlabel('Time (ms)')
set(gca, 'XGrid', 'on');

%saveas(gcf, ['../Analysis Output/EOG_subject_raster.svg'])
%close;


