
load('../mat_data/timelockeds/epochs_eog_all.mat');
load(['../mat_data/timelockeds/epochs_eog_avgrast.mat']);
load('../mat_data/timelockeds/epochs_eog.mat');
load('../mat_data/timelockeds/epochs_eog_resp.mat');


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


%% find biggest response from structure

epochs_eog_resp = struct;

%Calculate max response for subjects
%NB, min between 50-150ms (sample 111-131)
for i = 1:length(sub_date.ID)
    
    subinpath = ['../mat_data/timelockeds/ID' sub_date.ID{i} '/'];
    
    epochs_eog_resp.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:length(conditions)
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    label = eval(['cond.' char(conditions(ii)) 'label']);
    
        for stim_index = 1:nstim
        
        %window:
        minwin = cell2mat(epochs_eog.(conditions{ii})(i,stim_index));
        %minresponse
        minresp = min(minwin(111:131));
        
        epochs_eog_resp.(conditions{ii})(i, stim_index) = minresp
            
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

%To rectify vector a:
% a = [1 2 3 -4]
% a(a<0) = a(a<0)*-1

%[M, I] = min(mean(cell2mat(epochs_eog.PO60(:,6))));

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




%% Permutetest of EOG + manuscript plot

%create matrix (that works for permutest) of subject responses
for i = 1:numel(sub_date.ID)
    POavg(:,i) = epochs_eog.PO60{i,5}';
    GPavg(:,i) = epochs_eog.GP60{i,4}';
end

[clusters, p_values, t_sums, permutation_distribution] = permutest(POavg, GPavg, true, 0.01, 10000, true, inf);

%% Manuscript plots

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
plot(mean(POavg,2), 'Color', colors(5,:)); %Manually specify colors to match other plots
plot(mean(GPavg,2), 'Color', colors(10,:)); % -''-
plot(clusters{:}, repmat(0.5*10^-5, 1, length(clusters{:})), 'red', 'LineWidth', 3)

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

saveas(gcf, ['../Analysis Output/MEG_manus_EOG.svg']);
