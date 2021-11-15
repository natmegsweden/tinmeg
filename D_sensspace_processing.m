%% Load conditions for subjects, average trials via ft_timelockanalysis

%NB, fs_ds = 200 is hardcoded in for adapting TOI to gap in identifying amp and lat

all_cmb_avg = load('../mat_data/timelockeds/all_cmb_avg.mat');
all_cmb_avg = all_cmb_avg.all_cmb_avg;

gravg_cmb = load('../mat_data/timelockeds/grand_avg_cmb.mat');
gravg_cmb = gravg_cmb.gravg_cmb;

mean_sub = load(['../mat_data/timelockeds/mean_sub.mat']);
mean_sub = mean_sub.timelockeds_cmb;
mean_sub.time = round(mean_sub.time,3); %Round to avoid issues with floating point precision in plots

%top_chan are 6 unique highest response channels in combined grads in all PO60 conditions (i.e. unique(max_grad) except 'MEG2022+2023', 'MEG2042+2043' (superior parietal)
top_chan = {'MEG0242+0243', 'MEG1222+1223', 'MEG1322+1323', 'MEG1332+1333', 'MEG1342+1343', 'MEG1442+1443', 'MEG1522+1523', 'MEG1612+1613',  'MEG2422+2423', 'MEG2612+2613',  'MEG2642+2643'};

%% Process timelockeds per condition, gather subject- and grand average

all_cmb_avg = struct();
gravg_cmb = struct ();

for ii = 1%:length(conditions)

nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
trig = eval(['cond.' char(conditions(ii)) 'trig']);

    for iii = 5%1:nstim

        for i = 1:length(sub_date.ID)

        subinpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
        
        tempdat = load([subinpath conditions{ii} 'ica.mat']);
        tempdat = tempdat.([conditions{ii} 'ica']);
 

        cfg = [];
        cfg.covariance = 'no';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"
        cfg.preproc.demean = 'yes';
        
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
        
        %save individual timelockeds
        %save(['../mat_data/timelockeds/ID' sub_date.ID{i} '/' (cond.(([conditions{ii} 'label'])){iii}) '_tlks.mat'], 'timelockeds', '-v7.3');
        
        
        cfg = [];
%        timelockeds_cmb = ft_combineplanar(cfg, timelockeds);
        
%        all_cmb_avg.(conditions{ii}){i, iii} = timelockeds_cmb.avg;
        
        %Grab only 'avg' parameter in struct compatible with cat(mean())
        %all_avg.(conditions{ii}){i, iii} = timelockeds.avg;
        
        %Save some subject data structure for use with plot functions later
        %save('../mat_data/timelockeds/mean_sub.mat', 'timelockeds_cmb', '-v7.3');
        clear timelockeds timelockeds_cmb
        
        %Calculate grand average over trial and subject per condition
%        gravg_cmb.(conditions{ii}){iii} = mean(cat(3, all_cmb_avg.(conditions{ii}){:, iii}), 3);

        end

        %save timelockeds?
        
    end

end

clear i ii iii trig nstim

%save('../mat_data/timelockeds/grand_avg_cmb.mat', 'gravg_cmb', '-v7.3');
%save('../mat_data/timelockeds/all_cmb_avg.mat', 'all_cmb_avg', '-v7.3');

%% Butterfly with selected sensors

%Old, hand-picked sensors from MultiPlotER
%l_mag_chan = {'MEG1611', 'MEG1621', 'MEG1811', 'MEG1641', 'MEG1631', 'MEG1841', 'MEG1731', 'MEG1941', 'MEG1911'};
%r_mag_chan = {'MEG2421', 'MEG2411', 'MEG2221', 'MEG2431', 'MEG2441', 'MEG2231', 'MEG2511', 'MEG2321', 'MEG2311'};
%l_grad_chan = {'MEG1612+1613', 'MEG1622+1623', 'MEG1812+1813', 'MEG1642+1643', 'MEG1632+1633', 'MEG1842+1843', 'MEG1732+1733', 'MEG1942+1943', 'MEG1912+1913'};
%r_grad_chan = {'MEG2422+2423', 'MEG2412+2413', 'MEG2222+2223', 'MEG2432+2433', 'MEG2442+2443', 'MEG2232+2233', 'MEG2512+2513', 'MEG2322+2323', 'MEG2312+2313'};
%other_chan = {'MEG1332+1333', 'MEG1342+1343', 'MEG2612+2613', 'MEG0242+0243', 'MEG1322+1323'};

%GRAD-PLOTS PO60
for j = 1:6

cfg =[];
mean_sub.avg = gravg_cmb.PO60{j}

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
        if ismember(mean_sub.label{i}, top_chan)
        plot(mean_sub.avg(i,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);

        %Other GRADS in black (50% opacity)
        elseif mean_sub.label{i}(end) == int2str(3) %if GRAD - plot
        plot(mean_sub.avg(i,:), 'Color', [0 0 0 0.25]);
        
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
    
    %find max of first 102 chan (grads) at mean of samples 116:131 (75-150ms)
    [val, ind] = sort(mean(mean_sub.avg(1:102,116:131),2), 'descend');

    max_grad{j,1} = mean_sub.label{ind(1)};
    max_grad{j,2} = mean_sub.label{ind(2)};
    max_grad{j,3} = mean_sub.label{ind(3)};
    max_grad{j,4} = mean_sub.label{ind(4)};
    max_grad{j,5} = mean_sub.label{ind(5)};
    max_grad{j,6} = mean_sub.label{ind(6)};
    
    %saveas(gcf, ['../Analysis Output/' cond.PO60label{j} 'butterfly.svg']);
    %close
    
end

%GRAD-PLOTS GP60
for j = 1:4

cfg =[];
mean_sub.avg = gravg_cmb.GP60{j}

figure('Position', [400 400 1800 400]); hold on;

title(cond.GP60label{j}, 'Interpreter', 'none')

minylim = -0.5*10^-12;
maxylim = 10*10^-12;

xlim([13 165]); %length(mean_sub.time); [-440ms 320ms]
ylim([minylim maxylim]);

%TOI response
patch('Faces', [1 2 3 4], 'Vertices', [116 -0.5*10^-12; 116 10*10^-12; 131 10*10^-12; 131 -0.5*10^-12], 'FaceColor', 'green', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)

%An attempt to make defining shaded regions (baseline and gap) somewhat readable
if j == 1; b_low = find(mean_sub.time == -0.200); b_high = find(mean_sub.time == -0.0500);
            g_low = find(mean_sub.time == -0.0500); g_high = find(mean_sub.time == 0);
elseif j == 2; b_low = find(mean_sub.time == -0.260); b_high = find(mean_sub.time == -0.110);
            g_low = find(mean_sub.time == -0.110); g_high = find(mean_sub.time == -0.060);  
elseif j == 3; b_low = find(mean_sub.time == -0.320); b_high = find(mean_sub.time == -0.170);
            g_low = find(mean_sub.time == -0.170); g_high = find(mean_sub.time == -0.120);
elseif j == 4; b_low = find(mean_sub.time == -0.440); b_high = find(mean_sub.time == -0.290);
            g_low = find(mean_sub.time == -0.290); g_high = find(mean_sub.time == -0.240);
end

%Baseline shade
patch('Faces', [1 2 3 4], 'Vertices', [b_low minylim; b_low maxylim; b_high maxylim; b_high minylim], 'FaceColor', 'blue', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)
%Gap shade
patch('Faces', [1 2 3 4], 'Vertices', [g_low minylim; g_low maxylim; g_high maxylim; g_high minylim], 'FaceColor', 'red', 'FaceAlpha', 0.05, 'EdgeAlpha', 0)



    for i = 1:102
        
        %Top 6 channels (amplitude) in orange
        if ismember(mean_sub.label{i}, top_chan)
        plot(mean_sub.avg(i,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);

        %Other GRADS in black (50% opacity)
        elseif mean_sub.label{i}(end) == int2str(3) %if GRAD - plot
        plot(mean_sub.avg(i,:), 'Color', [0 0 0 0.25]);
        
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
    
    %saveas(gcf, ['../Analysis Output/' cond.GP60label{j} 'butterfly.svg']);
    %close
    
end

%GRAD-PLOTS GO60
for j = 1

cfg =[];
mean_sub.avg = gravg_cmb.GO{j}

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
        if ismember(mean_sub.label{i}, top_chan)
        plot(mean_sub.avg(i,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);

        %Other GRADS in black (50% opacity)
        elseif mean_sub.label{i}(end) == int2str(3) %if GRAD - plot
        plot(mean_sub.avg(i,:), 'Color', [0 0 0 0.25]);
        
        end
        
        %dashed line at time zero (sample 101)
        plot([101 101], [-0.5*10^-12 10*10^-12], 'k --');
        
        %Convert x-axis to ms from sample number
        x = gca;
        x.XTick = [1:10:165];
        x.XTickLabel = [-500:50:320];
        x.XMinorTick = 'on';
        x.FontSize = 20;
        x.XTickLabelRotation = 90;
    end
    
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

%% Quantify individual subjects components latencies and amplitude

%NB top_chan is defined during butterfly plots
top_chan_idx = find(ismember(mean_sub.label, top_chan));

sub_sensoi = struct();

%Gather mean of sensors of interest (from top_chan) for all subjects and trials
for i = 1:numel(conditions)
    n_stim = length(cond.([conditions{i} 'label']));
    
    for ii = 1:numel(sub_date.ID)
        
        for iii = 1:n_stim
           
        sub_sensoi.(conditions{i}){ii,iii} = mean(all_cmb_avg.(conditions{i}){ii,iii}(top_chan_idx,:));
            
        %for trial (column)
        end
    
    %for subject (row)
    end

%for condition
end

%save('../mat_data/timelockeds/subjects_sensoi_avg.mat', 'sub_sensoi', '-v7.3');

%% Inspect and plot PO60_90 (i.e. pulse only control). Also write latency of max and mean amplitude within TOI to struct.
%NB mean_sub loads up top

%Load average sensors of interest for subjects
sub_sensoi = load('../mat_data/timelockeds/subjects_sensoi_avg.mat');
sub_sensoi = sub_sensoi.sub_sensoi;

%Structure for amplitude and latency-measures
sub_amp_lat = struct();

%Specify plot y-limits
minylim = 0*10^-12;
maxylim = 14*10^-12;

%Window/patch colors
red = [241 88 84]/256;
orange = [250 164 58]/256;
green = [96 189 104]/256;
blue = [93 165 218]/256;
purple = [178 118 178]/256;

%number of subjects
n_subs = 1:numel(sub_date.ID);

%% PO60_90: Extract and plot amplitude and latencies

%Time windows of interest, varies with gap position!!
N1on = find(mean_sub.time == 0.050);
N1off = find(mean_sub.time == 0.150);
P2on = find(mean_sub.time == 0.150);
P2off = find(mean_sub.time == 0.250);

figure('Position', [400 400 1800 400]); hold on;
xlim([41 165]);
ylim([minylim maxylim]);
title('PO60_90', 'Interpreter', 'none');
for i = 1:numel(sub_date.ID)

%find lat and amp for N1, collect in struct
[M, I] = max(sub_sensoi.PO60{i,5}(N1on:N1off));
sub_amp_lat.PO60_90_N1lat(i,1) = mean_sub.time(N1on-1+I);
sub_amp_lat.PO60_90_N1amp(i,1) = mean(sub_sensoi.PO60{i,5}(N1on:N1off));
clear M I

%find lat and amp for P2, collect in struct
[M, I] = max(sub_sensoi.PO60{i,5}(P2on:P2off));
sub_amp_lat.PO60_90_P2lat(i,1) = mean_sub.time(P2on-1+I);
sub_amp_lat.PO60_90_P2amp(i,1) = mean(sub_sensoi.PO60{i,5}(P2on:P2off));
clear M I
    
plot(sub_sensoi.PO60{i,5}, 'Color', [0 0 0 0.5])
tempmean(i,1:165) = sub_sensoi.PO60{i,5};
end

%N1 patch
patch('Faces', [1 2 3 4], 'Vertices', [N1on minylim; N1on maxylim; N1off maxylim; N1off minylim], 'FaceColor', orange, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%P1 patch
patch('Faces', [1 2 3 4], 'Vertices', [P2on minylim; P2on maxylim; P2off maxylim; P2off minylim], 'FaceColor', green, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

plot(mean(tempmean), 'Color', [0.75 0 0], 'LineWidth', 1.5); clear tempmean;
plot([101 101], [minylim maxylim], 'k --');
x = gca;
x.XTick = [1:10:165];
x.XTickLabel = [-500:50:320];
x.XMinorTick = 'on';
x.FontSize = 20;
x.XTickLabelRotation = 90;

% saveas(gcf, ['../Analysis Output/PO60_90_subs.svg']);
% close;

%% GO_60: Extract and plot amplitude and latencies

%Time windows of interest, varies with gap position!!
N1on = find(mean_sub.time == 0.050);
N1off = find(mean_sub.time == 0.150);
P2on = find(mean_sub.time == 0.150);
P2off = find(mean_sub.time == 0.250);

gapon = find(mean_sub.time == 0);
gapoff = find(mean_sub.time == 0.05);

figure('Position', [400 400 1800 400]); hold on;
xlim([41 165]);
ylim([0*10^-12 14*10^-12]);
title('GO_60', 'Interpreter', 'none');
for i = 1:numel(sub_date.ID)
    
%find lat and amp for N1, collect in struct
[M, I] = max(sub_sensoi.GO{i,1}(N1on:N1off));
sub_amp_lat.GO60_90_N1lat(i,1) = mean_sub.time(N1on-1+I);
sub_amp_lat.GO60_90_N1amp(i,1) = mean(sub_sensoi.GO{i,1}(N1on:N1off));
clear M I

%find lat and amp for P2, collect in struct
[M, I] = max(sub_sensoi.GO{i,1}(P2on:P2off));
sub_amp_lat.GO60_90_P2lat(i,1) = mean_sub.time(P2on-1+I);
sub_amp_lat.GO60_90_P2amp(i,1) = mean(sub_sensoi.GO{i,1}(P2on:P2off));
clear M I
    
plot(sub_sensoi.GO{i,1}, 'Color', [0 0 0 0.5])
tempmean(i,1:165) = sub_sensoi.GO{i,1};
end

%Gap patch
patch('Faces', [1 2 3 4], 'Vertices', [gapon minylim; gapon maxylim; gapoff maxylim; gapoff minylim], 'FaceColor', red, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%N1 patch
patch('Faces', [1 2 3 4], 'Vertices', [N1on minylim; N1on maxylim; N1off maxylim; N1off minylim], 'FaceColor', orange, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%P2 patch
patch('Faces', [1 2 3 4], 'Vertices', [P2on minylim; P2on maxylim; P2off maxylim; P2off minylim], 'FaceColor', green, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

plot(mean(tempmean), 'Color', [0.75 0 0], 'LineWidth', 1.5); clear tempmean;
%plot([101 101], [minylim maxylim], 'k --');
x = gca;
x.XTick = [1:10:165];
x.XTickLabel = [-500:50:320];
x.XMinorTick = 'on';
x.FontSize = 20;
x.XTickLabelRotation = 90;

% saveas(gcf, ['../Analysis Output/GO60_subs.svg']);
% close;

%% GP60-trials: Extract and plot amplitude and latencies

%Pulse time window for i120 and i240
pulN1on = find(mean_sub.time == 0.050);
pulN1off = find(mean_sub.time == 0.150);

pulP2on = find(mean_sub.time == 0.150);
pulP2off = find(mean_sub.time == 0.250);

for ii = 1:numel(cond.GP60label) %order: i0, i60, i120, i240
    
    %Specify list of times needed to correct for adapting TOI to shifting ISI (i.e ISI + GAP duration)
    gapcomp = [0.050 0.110 0.170 0.290];

    %Adapt GAP N1 & P2 window to gap ISI
    N1on = find(mean_sub.time == round(0.050 - gapcomp(ii), 3));
    N1off = find(mean_sub.time == round(0.150 - gapcomp(ii), 3));
    P2on = find(mean_sub.time == round(0.150 - gapcomp(ii), 3));
    P2off = find(mean_sub.time == round(0.250 - gapcomp(ii), 3));

    gapon = find(mean_sub.time == -gapcomp(ii));
    gapoff = gapon + 10; %10 samples = 50ms, hardcoded fs_ds = 200;


    figure('Position', [400 400 1800 400]); hold on;
    xlim([41 165]);
    ylim([0*10^-12 14*10^-12]);
    title(cond.GP60label{ii}, 'Interpreter', 'none');
    for i = 1:numel(sub_date.ID)

    %find lat and amp for N1, collect in struct
    [M, I] = max(sub_sensoi.GP60{i,ii}(N1on:N1off));
    sub_amp_lat.([cond.GP60label{ii} '_N1lat'])(i,1) = mean_sub.time(N1on-1+I) + gapcomp(ii); %NB - Compensate t = 0 to first stimulation event (i.e gap onset)
    sub_amp_lat.([cond.GP60label{ii} '_N1amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(N1on:N1off));
    clear M I

    %find lat and amp for P2, collect in struct
    [M, I] = max(sub_sensoi.GP60{i,ii}(P2on:P2off));
    sub_amp_lat.([cond.GP60label{ii} '_P2lat'])(i,1) = mean_sub.time(P2on-1+I) + gapcomp(ii); %NB - Compensate t = 0 to first stimulation event (i.e gap onset)
    sub_amp_lat.([cond.GP60label{ii} '_P2amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(P2on:P2off));
    clear M I
        
        %For i120, pick out P2 response to Pulse
        if ii == 3;
            %find lat and amp for P2, collect in struct
            [M, I] = max(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off));
            sub_amp_lat.([cond.GP60label{ii} '_pulP2lat'])(i,1) = mean_sub.time(pulP2on-1+I);
            sub_amp_lat.([cond.GP60label{ii} '_pulP2amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off));
            clear M I
            
        elseif ii == 4;
            %find lat and amp for N1 and P2, collect in struct
            [M, I] = max(sub_sensoi.GP60{i,ii}(pulN1on:pulN1off));
            sub_amp_lat.([cond.GP60label{ii} '_pulN1lat'])(i,1) = mean_sub.time(pulN1on-1+I); %pulse is t = 0, no need to compensate
            sub_amp_lat.([cond.GP60label{ii} '_pulN1amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(pulN1on:pulN1off));
            clear M I
            
            [M, I] = max(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off));
            sub_amp_lat.([cond.GP60label{ii} '_pulP2lat'])(i,1) = mean_sub.time(pulP2on-1+I); %pulse is t = 0, no need to compensate
            sub_amp_lat.([cond.GP60label{ii} '_pulP2amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off));
            clear M I
            
        end

    plot(sub_sensoi.GP60{i,ii}, 'Color', [0 0 0 0.5])
    tempmean(i,1:165) = sub_sensoi.GP60{i,ii};
    end

    %Gap patch
    patch('Faces', [1 2 3 4], 'Vertices', [gapon minylim; gapon maxylim; gapoff maxylim; gapoff minylim], 'FaceColor', red, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    %N1 Patch
    patch('Faces', [1 2 3 4], 'Vertices', [N1on minylim; N1on maxylim; N1off maxylim; N1off minylim], 'FaceColor', orange, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    %P2 patch
    patch('Faces', [1 2 3 4], 'Vertices', [P2on minylim; P2on maxylim; P2off maxylim; P2off minylim], 'FaceColor', green, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
   
    if ii == 3; %i120
        %pulP2 patch
        patch('Faces', [1 2 3 4], 'Vertices', [pulP2on minylim; pulP2on maxylim; pulP2off maxylim; pulP2off minylim], 'FaceColor', purple, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
    elseif ii == 4; %i240
        %pulN1 patch
        patch('Faces', [1 2 3 4], 'Vertices', [pulN1on minylim; pulN1on maxylim; pulN1off maxylim; pulN1off minylim], 'FaceColor', blue, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);   
        %pulP2 patch
        patch('Faces', [1 2 3 4], 'Vertices', [pulP2on minylim; pulP2on maxylim; pulP2off maxylim; pulP2off minylim], 'FaceColor', purple, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);   
    end
    
    plot(mean(tempmean), 'Color', [0.75 0 0], 'LineWidth', 1.5); clear tempmean;
    plot([101 101], [minylim maxylim], 'k --');
    x = gca;
    x.XTick = [1:10:165];
    x.XTickLabel = [-500:50:320];
    x.XMinorTick = 'on';
    x.FontSize = 20;
    x.XTickLabelRotation = 90;

%     saveas(gcf, ['../Analysis Output/' cond.GP60label{ii} '_subs.svg']);
%     close;
    
    clear gapcomp

end

%% Plot amp and latencies

%save('../mat_data/timelockeds/sub_amp_lat.mat', 'sub_amp_lat', '-v7.3');

% sub_amp_lat = load('../mat_data/timelockeds/sub_amp_lat.mat');
% sub_amp_lat = sub_amp_lat.sub_amp_lat;

%Amplitude N1
figure('Position', [600 400 1200 600]);
subplot(1,2,1); hold on;
title('Amplitude N1', 'Interpreter', 'none');
xlim([0.75 7.25]);
ylim([0 9*10^-12]);

%N1
plot([1 2 3 4 5 6 7], [sub_amp_lat.PO60_90_N1amp sub_amp_lat.GO60_90_N1amp sub_amp_lat.GP60_i0_N1amp sub_amp_lat.GP60_i60_N1amp sub_amp_lat.GP60_i120_N1amp sub_amp_lat.GP60_i240_N1amp sub_amp_lat.GP60_i240_pulN1amp], 'Color', [0 0 0 0.5], 'HandleVisibility', 'off')
plot([1 2 3 4 5 6 7], [mean(sub_amp_lat.PO60_90_N1amp) mean(sub_amp_lat.GO60_90_N1amp) mean(sub_amp_lat.GP60_i0_N1amp) mean(sub_amp_lat.GP60_i60_N1amp) mean(sub_amp_lat.GP60_i120_N1amp) mean(sub_amp_lat.GP60_i240_N1amp) mean(sub_amp_lat.GP60_i240_pulN1amp)], 'Color', [0.75 0 0], 'LineWidth', 1.5, 'DisplayName', 'mean')

scatter(ones(1,22), sub_amp_lat.PO60_90_N1amp, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+1, sub_amp_lat.GO60_90_N1amp, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Gap N1')
scatter(ones(1,22)+2, sub_amp_lat.GP60_i0_N1amp, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+3, sub_amp_lat.GP60_i60_N1amp, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+4, sub_amp_lat.GP60_i120_N1amp, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+5, sub_amp_lat.GP60_i240_N1amp, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+6, sub_amp_lat.GP60_i240_pulN1amp, 'filled', 'MarkerFaceColor', blue, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Pulse N1')

ylabel('Mean amplitude within TOI');

x = gca;
x.XTick = [1:7]
x.XTickLabel = ({'Pulse Only', 'Gap Only', 'ISI 0ms', 'ISI 60ms', 'ISI120', 'ISI240', 'ISI240 Pulse'});
x.FontSize = 20;
x.XTickLabelRotation = -90;
x.YMinorTick = 'on';

legend show
legend boxoff

%Amplitude P2
subplot(1,2,2); hold on;
title('Amplitude P2', 'Interpreter', 'none');
xlim([0.75 8.25]);
ylim([0 9*10^-12]);

%P2
plot([1 2 3 4 5 6 7 8], [sub_amp_lat.PO60_90_P2amp sub_amp_lat.GO60_90_P2amp sub_amp_lat.GP60_i0_P2amp sub_amp_lat.GP60_i60_P2amp sub_amp_lat.GP60_i120_P2amp sub_amp_lat.GP60_i240_P2amp sub_amp_lat.GP60_i120_pulP2amp sub_amp_lat.GP60_i240_pulP2amp], 'Color', [0 0 0 0.5], 'HandleVisibility', 'off')
plot([1 2 3 4 5 6 7 8], [mean(sub_amp_lat.PO60_90_P2amp) mean(sub_amp_lat.GO60_90_P2amp) mean(sub_amp_lat.GP60_i0_P2amp) mean(sub_amp_lat.GP60_i60_P2amp) mean(sub_amp_lat.GP60_i120_P2amp) mean(sub_amp_lat.GP60_i240_P2amp) mean(sub_amp_lat.GP60_i120_pulP2amp) mean(sub_amp_lat.GP60_i240_pulP2amp)], 'Color', [0.75 0 0], 'LineWidth', 1.5, 'DisplayName', 'mean')

scatter(ones(1,22), sub_amp_lat.PO60_90_P2amp, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+1, sub_amp_lat.GO60_90_P2amp, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Gap P2')
scatter(ones(1,22)+2, sub_amp_lat.GP60_i0_P2amp, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+3, sub_amp_lat.GP60_i60_P2amp, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+4, sub_amp_lat.GP60_i120_P2amp, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+5, sub_amp_lat.GP60_i240_P2amp, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+6, sub_amp_lat.GP60_i120_pulP2amp, 'filled', 'MarkerFaceColor', purple, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Pulse P2')
scatter(ones(1,22)+7, sub_amp_lat.GP60_i240_pulP2amp, 'filled', 'MarkerFaceColor', purple, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')

x = gca;
x.XTick = [1:8]
x.XTickLabel = ({'Pulse Only', 'Gap Only', 'ISI 0ms', 'ISI 60ms', 'ISI120', 'ISI240', 'ISI120 Pulse', 'ISI240 Pulse'});
x.FontSize = 20;
x.XTickLabelRotation = -90;
x.YMinorTick = 'on';

legend show
legend boxoff

% saveas(gcf, ['../Analysis Output/amp_subs.svg']);
% close;

%Latency Plot
figure('Position', [600 400 1200 600]);
subplot(1,2,1); hold on;
title('Latency N1', 'Interpreter', 'none');
xlim([0.75 7.25]);
ylim([0 0.3]);

%N1
plot([1 2 3 4 5 6 7], [sub_amp_lat.PO60_90_N1lat sub_amp_lat.GO60_90_N1lat sub_amp_lat.GP60_i0_N1lat sub_amp_lat.GP60_i60_N1lat sub_amp_lat.GP60_i120_N1lat sub_amp_lat.GP60_i240_N1lat sub_amp_lat.GP60_i240_pulN1lat], 'Color', [0 0 0 0.5], 'HandleVisibility', 'off')
plot([1 2 3 4 5 6 7], [mean(sub_amp_lat.PO60_90_N1lat) mean(sub_amp_lat.GO60_90_N1lat) mean(sub_amp_lat.GP60_i0_N1lat) mean(sub_amp_lat.GP60_i60_N1lat) mean(sub_amp_lat.GP60_i120_N1lat) mean(sub_amp_lat.GP60_i240_N1lat) mean(sub_amp_lat.GP60_i240_pulN1lat)], 'Color', [0.75 0 0], 'LineWidth', 1.5, 'DisplayName', 'mean (n = 22)')

scatter(ones(1,22), sub_amp_lat.PO60_90_N1lat, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+1, sub_amp_lat.GO60_90_N1lat, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Gap N1')
scatter(ones(1,22)+2, sub_amp_lat.GP60_i0_N1lat, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+3, sub_amp_lat.GP60_i60_N1lat, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+4, sub_amp_lat.GP60_i120_N1lat, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+5, sub_amp_lat.GP60_i240_N1lat, 'filled', 'MarkerFaceColor', orange, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+6, sub_amp_lat.GP60_i240_pulN1lat, 'filled', 'MarkerFaceColor', blue, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Pulse N1')

ylabel('Latency (ms)');

x = gca;
x.XTick = [1:7]
x.XTickLabel = ({'Pulse Only', 'Gap Only', 'ISI 0ms', 'ISI 60ms', 'ISI120', 'ISI240', 'ISI240 Pulse'});
x.FontSize = 20;
x.XTickLabelRotation = -90;
x.YMinorTick = 'on';
x.YTickLabel = [0:50:300];

legend show
legend boxoff

%P2
subplot(1,2,2); hold on;
title('Latency P2', 'Interpreter', 'none');
xlim([0.75 8.25]);
ylim([0 0.3]);
plot([1 2 3 4 5 6 7 8], [sub_amp_lat.PO60_90_P2lat sub_amp_lat.GO60_90_P2lat sub_amp_lat.GP60_i0_P2lat sub_amp_lat.GP60_i60_P2lat sub_amp_lat.GP60_i120_P2lat sub_amp_lat.GP60_i240_P2lat sub_amp_lat.GP60_i120_pulP2lat sub_amp_lat.GP60_i240_pulP2lat], 'Color', [0 0 0 0.5], 'HandleVisibility', 'off')
plot([1 2 3 4 5 6 7 8], [mean(sub_amp_lat.PO60_90_P2lat) mean(sub_amp_lat.GO60_90_P2lat) mean(sub_amp_lat.GP60_i0_P2lat) mean(sub_amp_lat.GP60_i60_P2lat) mean(sub_amp_lat.GP60_i120_P2lat) mean(sub_amp_lat.GP60_i240_P2lat) mean(sub_amp_lat.GP60_i120_pulP2lat) mean(sub_amp_lat.GP60_i240_pulP2lat)], 'Color', [0.75 0 0], 'LineWidth', 1.5, 'DisplayName', 'mean')

scatter(ones(1,22), sub_amp_lat.PO60_90_P2lat, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+1, sub_amp_lat.GO60_90_P2lat, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Gap P2')
scatter(ones(1,22)+2, sub_amp_lat.GP60_i0_P2lat, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+3, sub_amp_lat.GP60_i60_P2lat, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+4, sub_amp_lat.GP60_i120_P2lat, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+5, sub_amp_lat.GP60_i240_P2lat, 'filled', 'MarkerFaceColor', green, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')
scatter(ones(1,22)+6, sub_amp_lat.GP60_i120_pulP2lat, 'filled', 'MarkerFaceColor', purple, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Pulse P2')
scatter(ones(1,22)+7, sub_amp_lat.GP60_i240_pulP2lat, 'filled', 'MarkerFaceColor', purple, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off')

x = gca;
x.XTick = [1:8]
x.XTickLabel = ({'Pulse Only', 'Gap Only', 'ISI 0ms', 'ISI 60ms', 'ISI120', 'ISI240', 'ISI120 Pulse', 'ISI240 Pulse'});
x.FontSize = 20;
x.XTickLabelRotation = -90;
x.YMinorTick = 'on';
x.YTickLabel = [0:50:300];

legend show
legend('Location', 'southeast');
legend boxoff

% saveas(gcf, ['../Analysis Output/lat_subs.svg']);
% close;

%% Export amps and lats as csv

sub_amp_lat_csv = struct2table(sub_amp_lat);

%writetable(sub_amp_lat_csv);


%% TopoplotER

%mean_sub, all_cmb_avg & gravg_cmb loads in butterfly plots

zlimlow = -15*10^-14;
zlimhigh = 15*10^-14;

%% GO N1
mean_sub.avg = gravg_cmb.GO{1,1};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

%if here!
cfg.xlim = [0.050 0.150];
cfg.baseline = [-0.100 0];

ft_topoplotER(cfg, mean_sub);
title('GO60 N1', 'Interpreter', 'none');

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

saveas(gcf, ['../Analysis Output/' 'topo_GO60N1.svg']);
close

%% GO P2
mean_sub.avg = gravg_cmb.GO{1,1};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

%if here!
cfg.xlim = [0.150 0.250];
cfg.baseline = [-0.100 0];

ft_topoplotER(cfg, mean_sub);
title('GO60 P2', 'Interpreter', 'none');

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

saveas(gcf, ['../Analysis Output/' 'topo_GO60P2.svg']);
close

%% PO60_90 N1
mean_sub.avg = gravg_cmb.PO60{1,5};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

%if here!
cfg.xlim = [0.050 0.150];
cfg.baseline = [-0.100 0];

ft_topoplotER(cfg, mean_sub);
title('PO60_90 N1', 'Interpreter', 'none');

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

saveas(gcf, ['../Analysis Output/' 'topo_PO60_90_N1.svg']);
close

%% PO60_90 P2
mean_sub.avg = gravg_cmb.PO60{5};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

%if here!
cfg.xlim = [0.150 0.250];
cfg.baseline = [-0.100 0];

ft_topoplotER(cfg, mean_sub);
title('PO60_90 P2', 'Interpreter', 'none');

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

saveas(gcf, ['../Analysis Output/' 'topo_PO60_90_P2.svg']);
close

%% Loop for GP60-trials

isicomp = [0 0.060 0.120 0.240]; %Shift x-lim so N1/P2 TOI match gap onset

ft_hastoolbox('brewermap', 1);   % ensure this toolbox is on the path

for i = 1:numel(cond.GP60label);

%N1
mean_sub.avg = gravg_cmb.GP60{1,i};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

cfg.xlim = [0-isicomp(i) 0.100-isicomp(i)];
cfg.baseline = [-0.150-isicomp(i) -0.050-isicomp(i)];

ft_topoplotER(cfg, mean_sub);
title(([cond.GP60label{i} ' N1']), 'Interpreter', 'none');

colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

saveas(gcf, ['../Analysis Output/topo_' cond.GP60label{i} '_N1.svg']);
close

%P2
cfg.xlim = [0.100-isicomp(i) 0.200-isicomp(i)];

ft_topoplotER(cfg, mean_sub);
title(([cond.GP60label{i} ' P2']), 'Interpreter', 'none');

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

% saveas(gcf, ['../Analysis Output/topo_' cond.GP60label{i} '_P2.svg']);
% close

end

%% Pulse components - WIP - Be careful - static!

i = 4;

%i240 N1
mean_sub.avg = gravg_cmb.GP60{1,i};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

cfg.xlim = [0.050 0.150];
cfg.baseline = [-0.390 -0.290];

ft_topoplotER(cfg, mean_sub);
title(([cond.GP60label{i} ' Pulse N1']), 'Interpreter', 'none');

colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

% saveas(gcf, ['../Analysis Output/topo_' cond.GP60label{i} '_pulN1.svg']);
% close

%i240 P2
mean_sub.avg = gravg_cmb.GP60{1,i};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

cfg.xlim = [0.150 0.250];
cfg.baseline = [-0.390 -0.290];

ft_topoplotER(cfg, mean_sub);
title(([cond.GP60label{i} ' Pulse P2']), 'Interpreter', 'none');

colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

% saveas(gcf, ['../Analysis Output/topo_' cond.GP60label{i} '_pulP2.svg']);
% close


i = 3;
%i120 P2
mean_sub.avg = gravg_cmb.GP60{1,i};

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

cfg.xlim = [0.150 0.250];
cfg.baseline = [-0.270 -0.170];

ft_topoplotER(cfg, mean_sub);
title(([cond.GP60label{i} ' Pulse P2']), 'Interpreter', 'none');

colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

% saveas(gcf, ['../Analysis Output/topo_' cond.GP60label{i} '_pulP2.svg']);
% close
