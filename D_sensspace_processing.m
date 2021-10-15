%% Load conditions for subjects, average trials via ft_timelockanalysis

all_avg = struct();
gravg = struct();

%% Process timelockeds per condition, gather subject- and grand average

for ii = 1%:length(conditions)

nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
trig = eval(['cond.' char(conditions(ii)) 'trig']);

    for iii = 1%:nstim

        for i = 1%:length(sub_date.ID)

        subinpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
        
        tempdat = load([subinpath conditions{ii} 'ica.mat']);
        tempdat = tempdat.([conditions{ii} 'ica']);
 

        cfg = [];
        cfg.covariance = 'no';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"
        cfg.preproc.demean = 'yes';
        
        %NB - Forgot GO below, re-run eventually
        
        if ismember(conditions{ii}, {'GO60', 'GO70', 'PO60', 'PO70'});
        %Baseline window: 150ms before pulse onset in PO trials
        cfg.preproc.baselinewindow = [-0.150 0];
        elseif ismember(conditions{ii}, {'GP60', 'GP70'});
            %Baseline window variable timepoint: 150ms before gap onset in GP trials
            if iii == 1 %ISI 0
            cfg.preproc.baselinewindow = [-0.200 -0.050];
            elseif iii == 2 %ISI 60
                cfg.preproc.baselinewindow = [-0.260 -0.110];
            elseif iii == 3 %ISI 120
                cfg.preproc.baselinewindow = [-0.320 -0.170];
            elseif iii == 4 %ISI 240
                cfg.preproc.baselinewindow = [-0.440 -0.290];
            end
        end
        
        cfg.preproc.lpfilter = 'yes';
        cfg.preproc.lpfreq = 70;
        cfg.trials = tempdat.trialinfo == trig(iii);
        
        timelockeds = ft_timelockanalysis(cfg, tempdat);
        
        clear tempdat

        %Grab only 'avg' parameter in struct compatible with cat(mean())
        all_avg.(conditions{ii}){i, iii} = timelockeds.avg;
        
        %save('../mat_data/timelockeds/mean_sub.mat', 'timelockeds', '-v7.3');
        clear timelockeds
        
        %Calculate grand average over trial and subject per condition
        gravg.(conditions{ii}){iii} = mean(cat(3, all_avg.(conditions{ii}){:, iii}), 3);

        end

        %save timelockeds?
        
    end

end

clear i ii iii trig nstim

%save('../mat_data/timelockeds/all_avg.mat', 'all_avg', '-v7.3');
%save('../mat_data/timelockeds/grand_avg.mat', 'gravg', '-v7.3');

%% Butterfly with selected sensors

load('../mat_data/timelockeds/all_avg.mat');
load('../mat_data/timelockeds/grand_avg.mat');

mean_sub = load(['../mat_data/timelockeds/mean_sub.mat']);
mean_sub = mean_sub.timelockeds; clear timelockeds;

%Old, hand-picked sensors from MultiPlotER
%l_mag_chan = {'MEG1611', 'MEG1621', 'MEG1811', 'MEG1641', 'MEG1631', 'MEG1841', 'MEG1731', 'MEG1941', 'MEG1911'};
%r_mag_chan = {'MEG2421', 'MEG2411', 'MEG2221', 'MEG2431', 'MEG2441', 'MEG2231', 'MEG2511', 'MEG2321', 'MEG2311'};
%l_grad_chan = {'MEG1612+1613', 'MEG1622+1623', 'MEG1812+1813', 'MEG1642+1643', 'MEG1632+1633', 'MEG1842+1843', 'MEG1732+1733', 'MEG1942+1943', 'MEG1912+1913'};
%r_grad_chan = {'MEG2422+2423', 'MEG2412+2413', 'MEG2222+2223', 'MEG2432+2433', 'MEG2442+2443', 'MEG2232+2233', 'MEG2512+2513', 'MEG2322+2323', 'MEG2312+2313'};
%other_chan = {'MEG1332+1333', 'MEG1342+1343', 'MEG2612+2613', 'MEG0242+0243', 'MEG1322+1323'};

%top_chan are 6 unique highest response channels in combined grads in all PO60 conditions (i.e. unique(max_grad) except 'MEG2022+2023', 'MEG2042+2043' (superior parietal)
top_chan = {'MEG0242+0243', 'MEG1222+1223', 'MEG1322+1323', 'MEG1332+1333', 'MEG1342+1343', 'MEG1442+1443', 'MEG1522+1523', 'MEG1612+1613',  'MEG2422+2423', 'MEG2612+2613',  'MEG2642+2643'};


%GRAD-PLOTS PO60 (Combines planar in loop)
for j = 1:6

cfg =[];
mean_sub.avg = gravg.PO60{j}

combined_sub = ft_combineplanar(cfg, mean_sub);

%Round to avoid issues with floating point precision in plots
combined_sub.time = round(combined_sub.time,3)

figure('Position', [400 400 1800 400]); hold on;

title(cond.PO60label{j}, 'Interpreter', 'none')
minylim = -0.5*10^-12;
maxylim = 10*10^-12;

xlim([13 165]); %length(mean_sub.time); [-440ms 320ms]
ylim([minylim maxylim]);

%TOI response
patch('Faces', [1 2 3 4], 'Vertices', [116 minylim; 116 maxylim; 131 maxylim; 131 minylim], 'FaceColor', 'green', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)

%Baseline
patch('Faces', [1 2 3 4], 'Vertices', [71 minylim; 71 maxylim; 101 maxylim; 101 minylim], 'FaceColor', 'blue', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)

    for i = 1:102
        
        %Top 6 channels (amplitude) in orange
        if ismember(combined_sub.label{i}, top_chan)
        plot(combined_sub.avg(i,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);

        %Other GRADS in black (50% opacity)
        elseif combined_sub.label{i}(end) == int2str(3) %if GRAD - plot
        plot(combined_sub.avg(i,:), 'Color', [0 0 0 0.25]);
        
        end
        
        %dashed line at time zero (sample 71)
        plot([101 101], [-0.5*10^-12 10*10^-12], 'k --');
        
        %Convert x-axis to ms from sample number
        x = gca;
        x.XTick = [1:10:165];
        x.XTickLabel = [-500:50:320];
        x.XMinorTick = 'on';
        x.FontSize = 20;
        x.XTickLabelRotation = 90;
    end
    
    %find max of first 102 chan (grads) at mean of samples 76:86 (25-75ms)
    [val, ind] = sort(mean(combined_sub.avg(1:102,116:131),2), 'descend');

    max_grad{j,1} = combined_sub.label{ind(1)};
    max_grad{j,2} = combined_sub.label{ind(2)};
    max_grad{j,3} = combined_sub.label{ind(3)};
    max_grad{j,4} = combined_sub.label{ind(4)};
    max_grad{j,5} = combined_sub.label{ind(5)};
    max_grad{j,6} = combined_sub.label{ind(6)};
    
    %saveas(gcf, ['../Analysis Output/' cond.PO60label{j} 'butterfly.svg']);
    %close
    
end

%GRAD-PLOTS GP60 (Combines planar in loop)
for j = 1:4

cfg =[];
mean_sub.avg = gravg.GP60{j}

combined_sub = ft_combineplanar(cfg, mean_sub);

%Round to avoid issues with floating point precision in plots
combined_sub.time = round(combined_sub.time,3)

figure('Position', [400 400 1800 400]); hold on;

title(cond.GP60label{j}, 'Interpreter', 'none')

minylim = -0.5*10^-12;
maxylim = 10*10^-12;

xlim([13 165]); %length(mean_sub.time); [-440ms 320ms]
ylim([minylim maxylim]);

%TOI response
patch('Faces', [1 2 3 4], 'Vertices', [116 -0.5*10^-12; 116 10*10^-12; 131 10*10^-12; 131 -0.5*10^-12], 'FaceColor', 'green', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)

%An attempt to make defining shaded regions (baseline and gap) somewhat readable
if j == 1; b_low = find(combined_sub.time == -0.200); b_high = find(combined_sub.time == -0.0500);
            g_low = find(combined_sub.time == -0.0500); g_high = find(combined_sub.time == 0);
elseif j == 2; b_low = find(combined_sub.time == -0.260); b_high = find(combined_sub.time == -0.110);
            g_low = find(combined_sub.time == -0.110); g_high = find(combined_sub.time == -0.060);  
elseif j == 3; b_low = find(combined_sub.time == -0.320); b_high = find(combined_sub.time == -0.170);
            g_low = find(combined_sub.time == -0.170); g_high = find(combined_sub.time == -0.120);
elseif j == 4; b_low = find(combined_sub.time == -0.440); b_high = find(combined_sub.time == -0.290);
            g_low = find(combined_sub.time == -0.290); g_high = find(combined_sub.time == -0.240);
end

%Baseline shade
patch('Faces', [1 2 3 4], 'Vertices', [b_low minylim; b_low maxylim; b_high maxylim; b_high minylim], 'FaceColor', 'blue', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)
%Gap shade
patch('Faces', [1 2 3 4], 'Vertices', [g_low minylim; g_low maxylim; g_high maxylim; g_high minylim], 'FaceColor', 'red', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)



    for i = 1:102
        
        %Top 6 channels (amplitude) in orange
        if ismember(combined_sub.label{i}, top_chan)
        plot(combined_sub.avg(i,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);

        %Other GRADS in black (50% opacity)
        elseif combined_sub.label{i}(end) == int2str(3) %if GRAD - plot
        plot(combined_sub.avg(i,:), 'Color', [0 0 0 0.25]);
        
        end
        
        %dashed line at time zero (sample 71)
        plot([101 101], [-0.5*10^-12 10*10^-12], 'k --');
        
        %Convert x-axis to ms from sample number
        x = gca;
        x.XTick = [1:10:165];
        x.XTickLabel = [-500:50:320];
        x.XMinorTick = 'on';
        x.FontSize = 20;
        x.XTickLabelRotation = 90;
    end
    
    %find max of first 102 chan (grads) at mean of samples 76:86 (25-75ms)
    [val, ind] = sort(mean(combined_sub.avg(1:102,116:131),2), 'descend');

    max_grad{j,1} = combined_sub.label{ind(1)};
    max_grad{j,2} = combined_sub.label{ind(2)};
    max_grad{j,3} = combined_sub.label{ind(3)};
    max_grad{j,4} = combined_sub.label{ind(4)};
    max_grad{j,5} = combined_sub.label{ind(5)};
    max_grad{j,6} = combined_sub.label{ind(6)};
    
    %saveas(gcf, ['../Analysis Output/' cond.GP60label{j} 'butterfly.svg']);
    %close
    
end

%GRAD-PLOTS GO60 (Combines planar in loop)
for j = 1

cfg =[];
mean_sub.avg = gravg.GO{j}

combined_sub = ft_combineplanar(cfg, mean_sub);

%Round to avoid issues with floating point precision in plots
combined_sub.time = round(combined_sub.time,3)

figure('Position', [400 400 1800 400]); hold on;

title(cond.GOlabel{j}, 'Interpreter', 'none')

minylim = -0.5*10^-12;
maxylim = 10*10^-12;

xlim([13 165]); %length(mean_sub.time); [-440ms 320ms]
ylim([minylim maxylim]);

%TOI response
%patch('Faces', [1 2 3 4], 'Vertices', [116 -0.5*10^-12; 116 10*10^-12; 131 10*10^-12; 131 -0.5*10^-12], 'FaceColor', 'green', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)

%Baseline shade
patch('Faces', [1 2 3 4], 'Vertices', [71 minylim; 71 maxylim; 101 maxylim; 101 minylim], 'FaceColor', 'blue', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)
%Gap shade
patch('Faces', [1 2 3 4], 'Vertices', [101 minylim; 101 maxylim; 111 maxylim; 111 minylim], 'FaceColor', 'red', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)


    for i = 1:102
        
        %Top 6 channels (amplitude) in orange
        if ismember(combined_sub.label{i}, top_chan)
        plot(combined_sub.avg(i,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);

        %Other GRADS in black (50% opacity)
        elseif combined_sub.label{i}(end) == int2str(3) %if GRAD - plot
        plot(combined_sub.avg(i,:), 'Color', [0 0 0 0.25]);
        
        end
        
        %dashed line at time zero (sample 71)
        plot([101 101], [-0.5*10^-12 10*10^-12], 'k --');
        
        %Convert x-axis to ms from sample number
        x = gca;
        x.XTick = [1:10:165];
        x.XTickLabel = [-500:50:320];
        x.XMinorTick = 'on';
        x.FontSize = 20;
        x.XTickLabelRotation = 90;
    end
    
    %find max of first 102 chan (grads) at mean of samples 76:86 (25-75ms)
    [val, ind] = sort(mean(combined_sub.avg(1:102,116:131),2), 'descend');

    max_grad{j,1} = combined_sub.label{ind(1)};
    max_grad{j,2} = combined_sub.label{ind(2)};
    max_grad{j,3} = combined_sub.label{ind(3)};
    max_grad{j,4} = combined_sub.label{ind(4)};
    max_grad{j,5} = combined_sub.label{ind(5)};
    max_grad{j,6} = combined_sub.label{ind(6)};
    
    %saveas(gcf, ['../Analysis Output/' cond.GOlabel{j} 'butterfly.svg']);
    %close
    
end


%% Sensshape plot with highlighted chips

cfg = [];
cfg.output = 'scalp';

mri_segmented = ft_volumesegment(cfg, template_mri);

%Transform to neuromag coordsys - same as for sensors
mri_segmented = ft_convert_coordsys(mri_segmented, 'neuromag');

%Create mesh skull
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);

%ft_plot_mesh(mesh_scalp);

%Plot higlighted sensors
colors = ones(306, 3);
%l_chips = {'MEG1612', 'MEG1622', 'MEG1812', 'MEG1642', 'MEG1632', 'MEG1842', 'MEG1732', 'MEG1942', 'MEG1912'};
%r_chips = {'MEG2422', 'MEG2412', 'MEG2222', 'MEG2432', 'MEG2442', 'MEG2232', 'MEG2512', 'MEG2322', 'MEG2312'};
%top_chips = {'MEG0242', 'MEG1322', 'MEG1332', 'MEG1342', 'MEG1442', 'MEG1612', 'MEG2422', 'MEG2522',  'MEG2612',  'MEG2642'}; %Other strong responses from GRADS (i.e. max_grad)
%l_chips_idx = ismember(mean_sub.label, l_chips);
%r_chips_idx = ismember(mean_sub.label, r_chips);
%top_chips_idx = ismember(mean_sub.label, top_chips);

%Index of the first gradiometer in top_chan before combine_planar to match plot
for j = 1:numel(top_chan)
top_chan_temp{j} = top_chan{j}(1:end-5);
end
top_chan_idx = ismember(mean_sub.label, top_chan_temp);
clear top_chan_temp

%Write colour vector to row in mean_sub.label (+1 to row for both grads on chip)
for i = 1:numel(mean_sub.label)
    
%     if l_chips_idx(i) == 1
%     colors(i, :) = [0 0 1];
%     colors(i+1, :) = [0 0 1];
%     elseif r_chips_idx(i) == 1
%     colors(i, :) = [1 0 0];
%     colors(i+1, :) = [1 0 0];
    if top_chan_idx(i) == 1
    colors(i, :) = [0.85 0.325 0.098];
    colors(i+1, :) = [0.85 0.325 0.098];
    end
    
end

%Move scalp to not clip through sensors and be positioned nicely in helmet
mesh_scalp.pos(:,3) = mesh_scalp.pos(:,3) - 10; %Down 10mm
mesh_scalp.pos(:,2) = mesh_scalp.pos(:,2) - 15; %Back 15mm

figure('Position', [400 200 1800 1000]);
hold on
subplot(1,2,1)
sensors = mean_sub.grad;
ft_plot_sens(sensors, 'facecolor', colors, 'facealpha', 0.7);

ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

subplot(1,2,2)
sensors = mean_sub.grad;
ft_plot_sens(sensors, 'facecolor', colors, 'facealpha', 0.7);

ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);

view([-100 25])

%% TopoplotER

load('../mat_data/sensspace/sensspace_avg.mat', 

for c = 1:numel(conditions);

mean_sub = sensspace_avg.(conditions{c}){1,1};

    cfg = [];
    cfg.parameter = 'avg';
    cfg.layout = 'neuromag306mag.lay';
    cfg.colorbar = 'no';
    cfg.zlim = [-7*10^-14 7*10^-14];
    
    %cfg.xlim = [0.002 0.3];
    %cfg.baseline = [-0.35 -0.1];
    
    cfg.comment = 'no';

    figure('Position', [300 500 1400 900]); %From: Left, Down - Width, Height
    ha = tight_subplot(1,6,[.05 .015],[.2 .1],[.05 .05]);
    
    for i = 1:numel(cond.([conditions{c} 'label'])); axes(ha(i))

        mean_sub.avg = gravg.(conditions{c}){i};

        %subplot(2,3,i);
        ft_topoplotER(cfg, mean_sub);
        %title(cond.PO60label{i}, 'Interpreter', 'none');

        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

    end
    
    saveas(gcf, ['../Analysis Output/' 'topoplot_' conditions{c} '.svg']);
    close
    
end
