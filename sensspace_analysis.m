%% Load conditions for subjects, average trials

%Load tlk_sub.mat -- see obsidian notes

%NB, fs_ds = 200 is hardcoded in for adapting TOI to gap in identifying amp and lat

all_cmb_avg = load('../mat_data/timelockeds/all_cmb_avg.mat');
all_cmb_avg = all_cmb_avg.all_cmb_avg;

gravg_cmb = load('../mat_data/timelockeds/grand_avg_cmb.mat');
gravg_cmb = gravg_cmb.gravg_cmb;

mean_sub = load(['../mat_data/timelockeds/mean_sub.mat']);
mean_sub = mean_sub.timelockeds_cmb;
mean_sub.time = round(mean_sub.time,3); %Round to avoid issues with floating point precision in plots

%Load sensor positions from one subjects data
sensors = load('../mat_data/preprocessing/ID0539/PO60_ds_clean.mat');
sensors = sensors.cleaned4mat.grad;

%find n = 8(?) top_chan in control condition (PO60/70 90)
%plot amplitude of topchans
%boxplot of top_chans
%helmet
%permutations test

%% Find 8 combined gradiometers with biggest response in grand average of subjects in PO60/70 90-trials

%How many combined grads of interest?
n_top_chan = 8;

%For 60dB Carrier
%Create empty cell array
max_grad = cell(1,n_top_chan);

%Overwrite mean_sub with grand average to maintain "labels"
mean_sub.avg = gravg_cmb.PO60{5}; %col 5 = PO6090

%find max of first 102 chan (grads) at mean of samples 116:131 (75-150ms)
[val, ind] = sort(mean(mean_sub.avg(1:102,116:131),2), 'descend');

%write top i cmb_grads to struct
for i = 1:n_top_chan
    max_grad{1,i} = mean_sub.label{ind(i)};
end

top_chan60 = unique(max_grad)';

%For 70dB Carrier
%Create empty cell array
max_grad = cell(1,n_top_chan);

%Overwrite mean_sub with grand average to maintain "labels"
mean_sub.avg = gravg_cmb.PO70{4}; %col 5 = PO7090

%find max of first 102 chan (grads) at mean of samples 116:131 (75-150ms)
[val, ind] = sort(mean(mean_sub.avg(1:102,116:131),2), 'descend');

%write top i cmb_grads to struct
for i = 1:n_top_chan
    max_grad{1,i} = mean_sub.label{ind(i)};
end

top_chan70 = unique(max_grad)';

%% Sensshape plot with highlighted chips

cfg = [];
cfg.output = 'scalp';

%Segment template MRI for head mesh
mri_segmented = ft_volumesegment(cfg, template_mri);

%Transform to neuromag coordsys - same as for sensors
mri_segmented = ft_convert_coordsys(mri_segmented, 'neuromag');

%Create mesh skull
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);

%Move scalp to not clip through sensors and be positioned nicely in helmet
mesh_scalp.pos(:,3) = mesh_scalp.pos(:,3) - 35; %Down
mesh_scalp.pos(:,2) = mesh_scalp.pos(:,2) - 15; %Back

%Plot higlighted sensors
%Write colour vector to row in mean_sub.label (+/-1 to row for all sensors on chip)
colors2 = ones(306, 3);

%Index based on the first gradiometer in top_chan
%Sensshape needs idx for all 306 sensors which complicates this a bit.
for i = 1:numel(top_chan60)
top_chan_temp{i} = top_chan60{i}(1:end-5);
end
top_chan_idx = ismember(sensors.label, top_chan_temp);
clear top_chan_temp

%match "colors" to index
for i = 1:numel(sensors.label)
    if top_chan_idx(i) == 1
    colors2(i, :) = [0.85 0.325 0.098];
    colors2(i+1, :) = [0.85 0.325 0.098];
    colors2(i-1, :) = [0.85 0.325 0.098];
    end
end

[top_chan60 top_chan70]
warning('Note that sensors are the same for PO60 and PO70 - 90dB when n of sensors = 8');

figure('Position', [400 200 1800 1000]); hold on
subplot(1,2,1)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

subplot(1,2,2)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([-100 25])

%% Sanity check figure of identified gradiometers

%Find index for gradiometers only
grad_idx = ismember(mean_sub.label, top_chan60);

figure; hold on;
for i = 1:102 %102 gradiometers
    if grad_idx(i) == 1
        plot(mean(gravg_cmb.PO60{5}(i,1:165),1), 'r');
    else plot(mean(gravg_cmb.PO60{5}(i,1:165),1), 'Color', [0 0 0 0.25]);
    end
end

%% Grand average for SOI only
gravg_soi = struct;

%Find index for gradiometers only
grad_idx = ismember(mean_sub.label, top_chan60);

%For each condition
for i = 1:numel(conditions);
    %Each stim
    for ii = 1:numel(cond.([conditions{i} 'label']))
        j=1; %keep track of rows in temp
        %Create temp struct of SOI
        for iii = 1:102 %102 gradiometers
            if grad_idx(iii) == 1
               temp(j,:) = gravg_cmb.(conditions{i}){ii}(iii,1:165);
               j=j+1;

               %Write mean of temp and clear
               gravg_soi.(conditions{i}){ii} = mean(temp,1);

            end
        end
        clear j temp
    end
end

%% Find biggest response per subject

sub_sens_amp = struct;

%Keep only soi per subject and find max of mean of soi
%NB, max between 50-150ms (sample 111-131)
for i = 1:numel(sub_date.ID)

    sub_sens_amp.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1:numel(conditions)
    
    nstim = numel(cond.([conditions{ii} 'trig']));
    
        for stim_index = 1:nstim
        
        temp = mean(all_cmb_avg.(conditions{ii}){i,stim_index}(grad_idx,:),1);
        maxresp = max(temp(111:131));

        sub_sens_amp.(conditions{ii})(i, stim_index) = maxresp;
        
        clear temp maxresp;
        %For stim
        end
    
    %For conditions
    end

%For subjects
end



%% Manuscript figures

% Inherited from raster plots
xtick = [1:20:165];
xticklab = [-500:100:320];
triglinex = [101 101];
xrange = [41 161];

txtsize = 10;

lineylims = [0 1.6*10^-11];
boxylims = [0 1.6*10^-11];

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

figure('Units', 'centimeters', 'Position',  [5 5 50 30]);
tiledlayout(3,4, 'TileSpacing','compact', 'Padding','compact');

%PO60
nexttile; hold on;
patch('Faces', [1 2 3 4], 'Vertices', [111 lineylims(1); 111 lineylims(2); 131 lineylims(2); 131 lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:6
    plot(gravg_soi.PO60{i}, 'Color', colors(i,:))
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

ylabel({"60dB carrier", "Gradiometer ERF"})

title({'Pulse only trials', 'Different pulse level'});

legend({"", "70", "75", "80", "85", "90", "95"}, 'Location', 'northwest'); %First one empty to skip patch
legend('boxoff');

nexttile; hold on;
boxplot(sub_sens_amp.PO60(:, :), [70:5:95], 'Symbol', 'ok');

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
patch('Faces', [1 2 3 4], 'Vertices', [111 lineylims(1); 111 lineylims(2); 131 lineylims(2); 131 lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:4
    plot(gravg_soi.GP60{i}, 'Color', colors(i+6,:))
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

legend({"", "0 ms", "60 ms", "120 ms", "240 ms"}, 'Location', 'northwest'); %First one empty to skip patch
legend('boxoff');

title({'Gap + Pulse trials', 'Different ISI'});

nexttile; boxplot(sub_sens_amp.GP60(:, :), [0 60 120 240], 'Symbol', 'ok');
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

%PO70
nexttile; hold on;
patch('Faces', [1 2 3 4], 'Vertices', [111 lineylims(1); 111 lineylims(2); 131 lineylims(2); 131 lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:5
    plot(gravg_soi.PO70{i}, 'Color', colors(i+1,:)) %+1 to color to match PO60
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
ylabel({"70dB carrier", "Gradiometer ERF"})
xlabel("Time (ms)");

legend({"", "75", "80", "85", "90", "95"}, 'Location', 'northwest'); %First one empty to skip patch
legend('boxoff');

%NaNs to pad missing 70dB pulse in 70dB carrier
nexttile; boxplot([repmat(NaN, 1, 22)' sub_sens_amp.PO70(:, 1:5)], [70:5:95], 'Symbol', 'ok');
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

%GP70
nexttile; hold on;
patch('Faces', [1 2 3 4], 'Vertices', [111 lineylims(1); 111 lineylims(2); 131 lineylims(2); 131 lineylims(1)], 'FaceColor', [0 0 0], 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
for i = 1:4
    plot(gravg_soi.GP70{i}, 'Color', colors(i+6,:))
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

legend({"", "0 ms", "60 ms", "120 ms", "240 ms"}, 'Location', 'northwest'); %First one empty to skip patch
legend('boxoff');

nexttile; boxplot(sub_sens_amp.GP70(:, :), [0 60 120 240], 'Symbol', 'ok');
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

nexttile;
nexttile;

nexttile; hold on;
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

nexttile; hold on;
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([-100 25])
