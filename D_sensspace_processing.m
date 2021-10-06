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
        
        %Different TOI for PO and GP compensating 50ms (gap length) for PO-trials
        %Use GP trials times-variable for plots/analysis
        if ismember(conditions{ii}, {'PO60', 'PO70'});
            toilow = -0.450;
            toihigh = 0.250;
        elseif ismember(conditions{ii}, {'GP60', 'GP70', 'GO'});
            toilow = -0.400;
            toihigh = 0.300;
        end;

        %Define toi for trials
        cfg = [];
        cfg.toilim = [toilow toihigh];

        tempdat = ft_redefinetrial(cfg, tempdat);

        cfg = [];
        cfg.covariance = 'no';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"
        cfg.preproc.demean = 'yes';
        
        %NB - Forgot GO below, re-run eventually
        
        if ismember(conditions{ii}, {'PO60', 'PO70'});
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
        
        %Legacy from earlier script
        %Collect all conditions averaged over trial for all subjects
        %sensspace_avg.(conditions{ii}){i, iii} = ft_timelockanalysis(cfg, tempdat);

        %Grab only 'avg' parameter in struct compatible with cat(mean())
        all_avg.(conditions{ii}){i, iii} = timelockeds.avg;
        
        save('../mat_data/sensspace/mean_sub.mat', 'timelockeds', '-v7.3');
        clear timelockeds
        
        %Calculate grand average over trial and subject per condition
        gravg.(conditions{ii}){iii} = mean(cat(3, all_avg.(conditions{ii}){:, iii}), 3);

        end

        %save timelockeds?
        
    end

end

%save('../mat_data/sensspace/all_avg.mat', 'all_avg', '-v7.3');
%save('../mat_data/sensspace/grand_avg.mat', 'gravg', '-v7.3');

%% Butterfly with selected sensors

load(['../mat_data/sensspace/grand_avg.mat']);
mean_sub = load(['../mat_data/sensspace/mean_sub.mat']);
mean_sub = mean_sub.timelockeds; clear timelockeds;

l_mag_chan = {'MEG1611', 'MEG1621', 'MEG1811', 'MEG1641', 'MEG1631', 'MEG1841', 'MEG1731', 'MEG1941', 'MEG1911'};
r_mag_chan = {'MEG2421', 'MEG2411', 'MEG2221', 'MEG2431', 'MEG2441', 'MEG2231', 'MEG2511', 'MEG2321', 'MEG2311'};

l_grad_chan = {'MEG1612+1613', 'MEG1622+1623', 'MEG1812+1813', 'MEG1642+1643', 'MEG1632+1633', 'MEG1842+1843', 'MEG1732+1733', 'MEG1942+1943', 'MEG1912+1913'};
r_grad_chan = {'MEG2422+2423', 'MEG2412+2413', 'MEG2222+2223', 'MEG2432+2433', 'MEG2442+2443', 'MEG2232+2233', 'MEG2512+2513', 'MEG2322+2323', 'MEG2312+2313'};

other_chan = {'MEG1332+1333', 'MEG1342+1343', 'MEG2612+2613', 'MEG0242+0243', 'MEG1322+1323'};

%top_chan are 6 unique highest response channels in combined grads in all PO60 conditions (i.e. unique(max_grad))
%old top_chan = {'MEG0242+0243', 'MEG1132+1133', 'MEG1322+1323', 'MEG1332+1333',  'MEG1342+1343', 'MEG1612+1613', 'MEG2422+2423', 'MEG2522+2523',  'MEG2612+2613', 'MEG2642+2643'};

top_chan = {'MEG0242+0243', 'MEG1322+1323', 'MEG1332+1333', 'MEG1342+1343', 'MEG1442+1443', 'MEG1612+1613', 'MEG2422+2423', 'MEG2522+2523',  'MEG2612+2613',  'MEG2642+2643'};

%GRAD-PLOTS (Combines planar in loop)
for j = 1%:6

cfg =[];
mean_sub.avg = gravg.PO60{j}

combined_sub = ft_combineplanar(cfg, mean_sub);

figure('Position', [400 400 1800 600]); hold on;

title(cond.PO60label{j}, 'Interpreter', 'none')
xlim([1 131]);
ylim([-0.5*10^-12 10*10^-12]);

    for i = 1:102

        %Left mag chips
%         if ismember(combined_sub.label{i}, l_grad_chan)
%         plot(combined_sub.avg(i,:), 'Color', [0 0 1], 'LineWidth', 1.5);
% 
%         %Right mag chips
%         elseif ismember(combined_sub.label{i}, r_grad_chan)
%         plot(combined_sub.avg(i,:), 'Color', [1 0 0], 'LineWidth', 1.5);
        
        %Top 6 channels (amplitude) in orange
        if ismember(combined_sub.label{i}, top_chan)
        plot(combined_sub.avg(i,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5);

        %Other GRADS in black (50% opacity)
        elseif combined_sub.label{i}(end) == int2str(3) %if GRAD - plot
        plot(combined_sub.avg(i,:), 'Color', [0 0 0 0.5]);
        
        end
        
        %dashed line at time zero (sample 71)
        %plot([71 71], [0 10*10^-12], 'k --');
        
        %Convert x-axis to ms from sample number
%         x = gca;
%         x.XTick = [1:5:131];
%         x.XTickLabel = [-350:25:300];
    end
    
    %find max of first 102 chan (grads) at mean of samples 76:86 (25-75ms)
    [val, ind] = sort(mean(combined_sub.avg(1:102,76:86),2), 'descend');

    max_grad{j,1} = combined_sub.label{ind(1)};
    max_grad{j,2} = combined_sub.label{ind(2)};
    max_grad{j,3} = combined_sub.label{ind(3)};
    max_grad{j,4} = combined_sub.label{ind(4)};
    max_grad{j,5} = combined_sub.label{ind(5)};
    max_grad{j,6} = combined_sub.label{ind(6)};
    
    %saveas(gcf, ['../Analysis Output/' cond.PO60label{j} 'butterfly.svg']);
    %close
    
end

%MAG-PLOTS not updated to match GRADS above
for j = 1:6

mean_sub.avg = gravg.PO60{j}
    
figure('Position', [400 400 1800 600]); hold on;

title(cond.PO60label{j}, 'Interpreter', 'none')
xlim([0 130]);
ylim([-3*10^-13 4*10^-13]);

    for i = 1:306

        %Left mag chips
        if mean_sub.label{i}(end) == int2str(1) & ismember(mean_sub.label{i}, l_mag_chan)
        plot(mean_sub.avg(i,:), 'Color', [0 0 1], 'LineWidth', 1.5);

        %Right mag chips
        elseif mean_sub.label{i}(end) == int2str(1) & ismember(mean_sub.label{i}, r_mag_chan)
        plot(mean_sub.avg(i,:), 'Color', [1 0 0], 'LineWidth', 1.5);

        %Other Mags
        elseif mean_sub.label{i}(end) == int2str(1) %if MAG - plot
        plot(mean_sub.avg(i,:), 'Color', [0 0 0 0.5]);
        
        end

        plot([71 71], [-3*10^-13 4*10^-13], 'k --');
        
    end

    %saveas(gcf, ['../Analysis Output/' cond.PO60label{j} 'butterfly.svg']);
    %close
    
end


%% GP60

l_mag_chan = {'MEG1611', 'MEG1621', 'MEG1811', 'MEG1641', 'MEG1631', 'MEG1841', 'MEG1731', 'MEG1941', 'MEG1911'};
r_mag_chan = {'MEG2421', 'MEG2411', 'MEG2221', 'MEG2431', 'MEG2441', 'MEG2231', 'MEG2511', 'MEG2321', 'MEG2311'};

for j = 1:4

mean_sub.avg = gravg.GP60{j}
    
figure('Position', [400 400 1800 600]); hold on;

title(cond.GP60label{j}, 'Interpreter', 'none')
xlim([0 130]);
ylim([-3*10^-13 4*10^-13]);

    for i = 1:306

        %Left mag chips
        if mean_sub.label{i}(end) == int2str(1) & ismember(mean_sub.label{i}, l_mag_chan)
        plot(mean_sub.avg(i,:), 'Color', [0 0 1], 'LineWidth', 1.5);

        %Right mag chips
        elseif mean_sub.label{i}(end) == int2str(1) & ismember(mean_sub.label{i}, r_mag_chan)
        plot(mean_sub.avg(i,:), 'Color', [1 0 0], 'LineWidth', 1.5);

        %Other Mags
        elseif mean_sub.label{i}(end) == int2str(1) %if MAG - plot
        plot(mean_sub.avg(i,:), 'Color', [0 0 0 0.5]);
        
        end
        
        plot([71 71], [-3*10^-13 4*10^-13], 'k --');
        
    end
    
    saveas(gcf, ['../Analysis Output/' cond.GP60label{j} 'butterfly.svg']);
    close
    
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
cfg.numvertices = 800;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);

%ft_plot_mesh(mesh_scalp);

%Plot higlighted sensors
colors = ones(306, 3);
l_chips = {'MEG1612', 'MEG1622', 'MEG1812', 'MEG1642', 'MEG1632', 'MEG1842', 'MEG1732', 'MEG1942', 'MEG1912'};
r_chips = {'MEG2422', 'MEG2412', 'MEG2222', 'MEG2432', 'MEG2442', 'MEG2232', 'MEG2512', 'MEG2322', 'MEG2312'};

top_chips = {'MEG0242', 'MEG1322', 'MEG1332', 'MEG1342', 'MEG1442', 'MEG1612', 'MEG2422', 'MEG2522',  'MEG2612',  'MEG2642'}; %Other strong responses from GRADS (i.e. max_grad)

l_chips_idx = ismember(mean_sub.label, l_chips);
r_chips_idx = ismember(mean_sub.label, r_chips);

top_chips_idx = ismember(mean_sub.label, top_chips);

%Write colour vector to row in mean_sub.label (+1 to row for both grads on chip)
for i = 1:numel(mean_sub.label)
    
%     if l_chips_idx(i) == 1
%     colors(i, :) = [0 0 1];
%     colors(i+1, :) = [0 0 1];
%     elseif r_chips_idx(i) == 1
%     colors(i, :) = [1 0 0];
%     colors(i+1, :) = [1 0 0];
    if top_chips_idx(i) == 1
    colors(i, :) = [0.85 0.325 0.098];
    colors(i+1, :) = [0.85 0.325 0.098];
    end
    
end

figure('Position', [400 200 1800 1000]);
hold on
subplot(1,2,1)
sensors = mean_sub.grad;
ft_plot_sens(sensors, 'facecolor', colors, 'facealpha', 0.8);

ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [0.5 0.5 0.5]);
view([100 25])

subplot(1,2,2)
sensors = mean_sub.grad;
ft_plot_sens(sensors, 'facecolor', colors, 'facealpha', 0.8);

ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [0.5 0.5 0.5]);

view([-100 25])

%% MultiplotER

%load('../mat_data/sensspace/grand_avg.mat');

%Put grand average in a subjects structure for plotting
mean_sub = PO60_90{1,1};
mean_sub.avg = gravg.PO60{5};

%GRADS
cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306planar.lay';
cfg.showcomment = 'no';
cfg.channel = 'meggrad';
cfg.linecolor = [repmat(.5, 22, 3); 1 0 0]; %change according to n of subjects

ft_multiplotER(cfg, PO60_90{1:22}, mean_sub);

%MAGS
cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.showcomment = 'no';
cfg.channel = 'megmag';
cfg.linecolor = [repmat(.5, 22, 3); 1 0 0]; %change according to n of subjects

ft_multiplotER(cfg, PO60_90{1:22}, mean_sub);

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
