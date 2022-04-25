
txtsize = 16;

%% Load and plot relevant EOG

load(['../mat_data/timelockeds/epochs_eog_avgrast.mat']);

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
ha = tight_subplot(2,3,[.05 .015],[.2 .1],[.05 .05]);
colbar = [-8*10^-5 8*10^-5]; %set limits of color-gradient/"z-axis"

axes(ha(3));
colorbar('westoutside', 'AxisLocation', 'in');
ax = gca;
ax.FontSize = txtsize;
ax.CLim = colbar;

% 1 - PO60
    axes(ha(1));
    imagesc(epochs_eog_avgrast.PO60{1,5}, colbar);
    colormap parula;
    hold on

    plot(triglinex, [0 45], 'k --');

    ax = gca;

    ax.YTick = [ytick];
    ax.YTickLabel = [yticklab];
    %ax.YGrid = 'On';

    ax.XGrid = 'On';

    xlim(xrange);

    ax.XTick = [xtick];
    ax.XTickLabel = [];

    ax.FontSize = txtsize;
    
% 2 - G+P
    axes(ha(2));
    imagesc(epochs_eog_avgrast.GP60{1,4}, colbar);
    colormap parula;
    hold on

    plot(triglinex, [0 45], 'k --');

    ax = gca;

    ax.YTick = [ytick];
    ax.YTickLabel = [];
    %ax.YGrid = 'On';

    ax.XGrid = 'On';

    xlim(xrange);

    ax.XTick = [xtick];
    ax.XTickLabel = [];

    ax.FontSize = txtsize;


% 3 - Grand Average Amp PO
    axes(ha(4));
    plot(mean(epochs_eog_avgrast.PO60{1,5}), 'k', 'LineWidth', 0.5);
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

    ax.XTickLabel(2,:) = nan;
    ax.XTickLabel(4,:) = nan;
    ax.XTickLabel(6,:) = nan;
    ax.XTickLabel(8,:) = nan;
    ax.XTickLabel(10,:) = nan;
    ax.XTickLabel(12,:) = nan;
    ax.XTickLabel(14,:) = nan;
    ax.XTickLabel(16,:) = nan;

    ax.FontSize = txtsize;
    
% 4 - Grand Average Amp G+P
    axes(ha(5));
    plot(mean(epochs_eog_avgrast.GP60{1,4}), 'k', 'LineWidth', 0.5);
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

    ax.XTickLabel(2,:) = nan;
    ax.XTickLabel(4,:) = nan;
    ax.XTickLabel(6,:) = nan;
    ax.XTickLabel(8,:) = nan;
    ax.XTickLabel(10,:) = nan;
    ax.XTickLabel(12,:) = nan;
    ax.XTickLabel(14,:) = nan;
    ax.XTickLabel(16,:) = nan;
    
    ax.YTickLabel = [];

    ax.FontSize = txtsize;
    
saveas(gcf, ['../Analysis Output/Poster_EOG.svg']);
close;
    
%% Sensor space plots

mean_sub = load(['../mat_data/timelockeds/mean_sub.mat']);
mean_sub = mean_sub.timelockeds_cmb;
mean_sub.time = round(mean_sub.time,3); %Round to avoid issues with floating point precision in plots

%top_chan are 6 unique highest response channels in combined grads in all PO60 conditions (i.e. unique(max_grad) except 'MEG2022+2023', 'MEG2042+2043' (superior parietal)
top_chan = {'MEG0242+0243', 'MEG1222+1223', 'MEG1322+1323', 'MEG1332+1333', 'MEG1342+1343', 'MEG1442+1443', 'MEG1522+1523', 'MEG1612+1613',  'MEG2422+2423', 'MEG2612+2613',  'MEG2642+2643'};

%Load average sensors of interest for subjects
sub_sensoi = load('../mat_data/timelockeds/subjects_sensoi_avg.mat');
sub_sensoi = sub_sensoi.sub_sensoi;


%% Plots

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
n_subs = 22;

%PO60_90

    %Time windows of interest, varies with gap position!!
    N1on = find(mean_sub.time == 0.050);
    N1off = find(mean_sub.time == 0.150);
    P2on = find(mean_sub.time == 0.150);
    P2off = find(mean_sub.time == 0.250);

    figure('Position', [400 400 1800 400]); hold on;
    xlim([41 165]);
    ylim([minylim maxylim]);
    title('PO60_90', 'Interpreter', 'none');
    
    for i = 1:n_subs

    %find lat and amp for N1, collect in struct
    [M, I] = max(sub_sensoi.PO60{i,5}(N1on:N1off));
    sub_amp_lat.PO60_90_N1lat(i,1) = mean_sub.time(N1on-1+I);
    sub_amp_lat.PO60_90_N1amp(i,1) = mean(sub_sensoi.PO60{i,5}(N1on:N1off)); %Mean amp
    sub_amp_lat.PO60_90_N1amp_peak(i,1) = max(sub_sensoi.PO60{i,5}(N1on:N1off)); %Peak amp
    clear M I

    %find lat and amp for P2, collect in struct
    [M, I] = max(sub_sensoi.PO60{i,5}(P2on:P2off));
    sub_amp_lat.PO60_90_P2lat(i,1) = mean_sub.time(P2on-1+I);
    sub_amp_lat.PO60_90_P2amp(i,1) = mean(sub_sensoi.PO60{i,5}(P2on:P2off)); %Mean amp
    sub_amp_lat.PO60_90_P2amp_peak(i,1) = max(sub_sensoi.PO60{i,5}(P2on:P2off)); %Peak amp
    clear M I

    plot(sub_sensoi.PO60{i,5}, 'Color', [0 0 0 0.5])
    tempmean(i,1:165) = sub_sensoi.PO60{i,5};
    end

    %N1 patch
    patch('Faces', [1 2 3 4], 'Vertices', [N1on minylim; N1on maxylim; N1off maxylim; N1off minylim], 'FaceColor', blue, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    %P1 patch
    patch('Faces', [1 2 3 4], 'Vertices', [P2on minylim; P2on maxylim; P2off maxylim; P2off minylim], 'FaceColor', purple, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    plot(mean(tempmean), 'Color', [0.75 0 0], 'LineWidth', 1.5); %clear tempmean;
    plot([101 101], [minylim maxylim], 'k --');
    x = gca;
    x.XTick = [1:10:165];
    x.XTickLabel = [-500:50:320];
    x.XMinorTick = 'on';
    x.FontSize = txtsize;
    x.XTickLabelRotation = 90;

    saveas(gcf, ['../Analysis Output/Poster_PO60.svg']);
    close;
    
%GP

    %Pulse time window for i120 and i240
    pulN1on = find(mean_sub.time == 0.050);
    pulN1off = find(mean_sub.time == 0.150);

    pulP2on = find(mean_sub.time == 0.150);
    pulP2off = find(mean_sub.time == 0.250);

    for ii = 4

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
        for i = 1:n_subs;

        %find lat and amp for N1, collect in struct
        [M, I] = max(sub_sensoi.GP60{i,ii}(N1on:N1off));
        sub_amp_lat.([cond.GP60label{ii} '_N1lat'])(i,1) = mean_sub.time(N1on-1+I) + gapcomp(ii); %NB - Compensate t = 0 to first stimulation event (i.e gap onset)
        sub_amp_lat.([cond.GP60label{ii} '_N1amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(N1on:N1off)); %mean amp
        sub_amp_lat.([cond.GP60label{ii} '_N1amp_peak'])(i,1) = max(sub_sensoi.GP60{i,ii}(N1on:N1off)); %peak amp
        clear M I

        %find lat and amp for P2, collect in struct
        [M, I] = max(sub_sensoi.GP60{i,ii}(P2on:P2off));
        sub_amp_lat.([cond.GP60label{ii} '_P2lat'])(i,1) = mean_sub.time(P2on-1+I) + gapcomp(ii); %NB - Compensate t = 0 to first stimulation event (i.e gap onset)
        sub_amp_lat.([cond.GP60label{ii} '_P2amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(P2on:P2off)); %mean amp
        sub_amp_lat.([cond.GP60label{ii} '_P2amp_peak'])(i,1) = max(sub_sensoi.GP60{i,ii}(P2on:P2off)); %peak amp
        clear M I

            %For i120, pick out P2 response to Pulse
            if ii == 3;
                %find lat and amp for P2, collect in struct
                [M, I] = max(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off));
                sub_amp_lat.([cond.GP60label{ii} '_pulP2lat'])(i,1) = mean_sub.time(pulP2on-1+I);
                sub_amp_lat.([cond.GP60label{ii} '_pulP2amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off)); %mean amp
                sub_amp_lat.([cond.GP60label{ii} '_pulP2amp_peak'])(i,1) = max(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off)); %peak amp
                clear M I

            elseif ii == 4;
                %find lat and amp for N1 and P2, collect in struct
                [M, I] = max(sub_sensoi.GP60{i,ii}(pulN1on:pulN1off));
                sub_amp_lat.([cond.GP60label{ii} '_pulN1lat'])(i,1) = mean_sub.time(pulN1on-1+I); %pulse is t = 0, no need to compensate
                sub_amp_lat.([cond.GP60label{ii} '_pulN1amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(pulN1on:pulN1off)); %mean amp
                sub_amp_lat.([cond.GP60label{ii} '_pulN1amp_peak'])(i,1) = max(sub_sensoi.GP60{i,ii}(pulN1on:pulN1off)); %peak amp
                clear M I

                [M, I] = max(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off));
                sub_amp_lat.([cond.GP60label{ii} '_pulP2lat'])(i,1) = mean_sub.time(pulP2on-1+I); %pulse is t = 0, no need to compensate
                sub_amp_lat.([cond.GP60label{ii} '_pulP2amp'])(i,1) = mean(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off)); %mean amp
                sub_amp_lat.([cond.GP60label{ii} '_pulP2amp_peak'])(i,1) = max(sub_sensoi.GP60{i,ii}(pulP2on:pulP2off)); %peak amp
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

        plot(mean(tempmean), 'Color', [0.75 0 0], 'LineWidth', 1.5); %clear tempmean;
        plot([101 101], [minylim maxylim], 'k --');
        x = gca;
        x.XTick = [1:10:165];
        x.XTickLabel = [-500:50:320];
        x.XMinorTick = 'on';
        x.FontSize = txtsize;
        x.XTickLabelRotation = 90;

        saveas(gcf, ['../Analysis Output/Poster_GP_subs.svg']);
        close;

        clear gapcomp

    end
 
 %% Dips
 
%Standard headmodel from running standard_mri throught MRI_prep2 pipeline
load('../mat_data/MRI_mat/standard_headmodel.mat');

headmodel_std = ft_convert_units(headmodel_std, 'cm');

%NB - Load dip positions - 3 in filename is updated
all_dip_pos = load('../mat_data/source_reconstruction/all_dip_pos3.mat');
all_dip_pos = all_dip_pos.all_dip_pos;

%% Mesh brain
figure('Renderer', 'painters', 'Position', [400 400 1800 1800]); hold on;
ft_plot_headmodel(headmodel_std, 'facealpha', 0.2, 'edgecolor', [0.5 0.5 0.5]);

for j = 1:22
   
    if all_dip_pos.PO60_90Lsrc(j,1) < 0
    hemicol = [0 0.75 1];
    elseif all_dip_pos.PO60_90Lsrc(j,1) > 0
    hemicol = [1 0 0];
    end
    
    ft_plot_dipole(all_dip_pos.PO60_90Lsrc(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', hemicol, 'alpha', 0.5); hold on
    
       
    if all_dip_pos.PO60_90Rsrc(j,1) < 0
    hemicol = [0 0.75 1];
    elseif all_dip_pos.PO60_90Rsrc(j,1) > 0
    hemicol = [1 0 0];
    end
    
    ft_plot_dipole(all_dip_pos.PO60_90Rsrc(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', hemicol, 'alpha', 0.5); hold on  
    
end
view([0.5 0.75 0])

saveas(gcf, ['../Analysis Output/Poster_dippos.svg']);
close;

%% Virt chans

all_virtch = load(['../mat_data/source_reconstruction/all_virtch_0-70.mat']);
all_virtch = all_virtch.all_virtch;


%Hacky time vector for Fs = 200;
timevec = [-0.500:0.005:0.320];
timevec = round(timevec, 3);

%Move virtual channels to the correct hemisphere, i.e. so that x coordinate > 0 is in right hemisphere and x < 0 is left
coi = {'PO60_90', 'GO_60', 'GP60_i0', 'GP60_i60', 'GP60_i120', 'GP60_i240'};

for ii = 1:numel(coi);

    for i = 1:n_subs

        %Only using PO as source position for virtual channel so only check this trial
        if all_dip_pos.PO60_90Rsrc(i,1) < 0;

           tempL = all_virtch.([coi{ii} 'R'])(i,:);
           all_virtch.([coi{ii} 'R'])(i,:) = all_virtch.([coi{ii} 'L'])(i,:);
           all_virtch.([coi{ii} 'L'])(i,:) = tempL;

        end

    end

end

clear coi

%Window/patch colors
red = [241 88 84]/256;
orange = [250 164 58]/256;
green = [96 189 104]/256;
blue = [93 165 218]/256;
purple = [178 118 178]/256;

minylim = -2.5*10^-11;
maxylim = 1.5*10^-11;

minxlim = 41;
maxxlim = 165;

%% 
figure('Position', [400 400 1600 800]); hold on;

% Pulse only
    %SUBPLOT R
    subplot(2,1,1); hold on;
    %title('Pulse, Right virtual channel response', 'Interpreter', 'none');
    for i = 1:n_subs;
        plot(all_virtch.PO60_90R(i,:), 'Color', [0 0 0 0.5])
        tempmeanR(i,1:165) = all_virtch.PO60_90R(i,:);
    end

    plot(mean(tempmeanR), 'Color', [0.75 0 0], 'LineWidth', 1.5); clear tempmean;
    plot([101 101], [minylim maxylim], 'k --');

    %N1 patch
    patch('Faces', [1 2 3 4], 'Vertices', [111 minylim; 111 maxylim; 131 maxylim; 131 minylim], 'FaceColor', blue, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    %P2 patch
    patch('Faces', [1 2 3 4], 'Vertices', [131 minylim; 131 maxylim; 151 maxylim; 151 minylim], 'FaceColor', purple, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    %Formatting
    xlim([minxlim maxxlim]);
    ylim([minylim maxylim]);

    x = gca;
    x.XTick = [1:10:165];
    x.XTickLabel = [];
    x.XMinorTick = 'on';
    x.FontSize = txtsize;
    x.XTickLabelRotation = 90;


    %SUBPLOT L
    subplot(2,1,2); hold on;
    %title('Pulse, Left virtual channel response', 'Interpreter', 'none');
    for i = 1:n_subs
    plot(all_virtch.PO60_90L(i,:), 'Color', [0 0 0 0.5])
    tempmeanL(i,1:165) = all_virtch.PO60_90L(i,:);
    end
    plot(mean(tempmeanL), 'Color', [0 0 0.75], 'LineWidth', 1.5); clear tempmean;
    plot([101 101], [minylim maxylim], 'k --');

    %N1 patch
    patch('Faces', [1 2 3 4], 'Vertices', [111 minylim; 111 maxylim; 131 maxylim; 131 minylim], 'FaceColor', blue, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    %P2 patch
    patch('Faces', [1 2 3 4], 'Vertices', [131 minylim; 131 maxylim; 151 maxylim; 151 minylim], 'FaceColor', purple, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

    %Formatting
    xlim([minxlim maxxlim]);
    ylim([minylim maxylim]);

    x = gca;
    x.XTick = [1:10:165];
    x.XTickLabel = [-500:50:320];
    x.XMinorTick = 'on';
    x.FontSize = txtsize;
    x.XTickLabelRotation = 90;

    saveas(gcf, ['../Analysis Output/Poster_virtchan_PO.svg']);
    close;
    
%% Gap + Pulse

%Pulse time window for i120 and i240
pulN1on = find(timevec == 0.050);
pulN1off = find(timevec == 0.150);

pulP2on = find(timevec == 0.150);
pulP2off = find(timevec == 0.250);

for ii = 4

    %Specify list of times needed to correct for adapting TOI to shifting ISI (i.e ISI + GAP duration)
    gapcomp = [0.050 0.110 0.170 0.290];

    %Adapt GAP N1 & P2 window to gap ISI
    N1on = find(timevec == round(0.050 - gapcomp(ii), 3));
    N1off = find(timevec == round(0.150 - gapcomp(ii), 3));
    P2on = find(timevec == round(0.150 - gapcomp(ii), 3));
    P2off = find(timevec == round(0.250 - gapcomp(ii), 3));

    gapon = find(timevec == -gapcomp(ii));
    gapoff = gapon + 10; %10 samples = 50ms, hardcoded fs_ds = 200;
    
    figure('Position', [400 400 1600 800]); hold on;

        %SUBPLOT R
        subplot(2,1,1); hold on;
        %title([cond.GP60label{ii} ', Right virtual channel response'], 'Interpreter', 'none');
        for i = 1:n_subs
            plot(all_virtch.([cond.GP60label{ii} 'R'])(i,:), 'Color', [0 0 0 0.5])
            tempmeanR(i,1:165) = all_virtch.([cond.GP60label{ii} 'R'])(i,:);
        end

        plot(mean(tempmeanR), 'Color', [0.75 0 0], 'LineWidth', 1.5); clear tempmean;
        plot([101 101], [minylim maxylim], 'k --');

        %Gap patch
        patch('Faces', [1 2 3 4], 'Vertices', [gapon minylim; gapon maxylim; gapoff maxylim; gapoff minylim], 'FaceColor', red, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
        %N1 patch
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
        
        %Formatting
        xlim([minxlim maxxlim]);
        ylim([minylim maxylim]);

        x = gca;
        x.XTick = [1:10:165];
        x.XTickLabel = [];
        x.XMinorTick = 'on';
        x.FontSize = txtsize;
        x.XTickLabelRotation = 90;


        %SUBPLOT L
        subplot(2,1,2); hold on;
        %title([cond.GP60label{ii} ', Left virtual channel response'], 'Interpreter', 'none');
        for i = 1:n_subs
        plot(all_virtch.([cond.GP60label{ii} 'L'])(i,:), 'Color', [0 0 0 0.5])
        tempmeanL(i,1:165) = all_virtch.([cond.GP60label{ii} 'L'])(i,:);
        end
        plot(mean(tempmeanL), 'Color', [0 0 0.75], 'LineWidth', 1.5); clear tempmean;
        plot([101 101], [minylim maxylim], 'k --');

        %Gap patch
        patch('Faces', [1 2 3 4], 'Vertices', [gapon minylim; gapon maxylim; gapoff maxylim; gapoff minylim], 'FaceColor', red, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
        %N1 patch
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

        %Formatting
        xlim([minxlim maxxlim]);
        ylim([minylim maxylim]);

        x = gca;
        x.XTick = [1:10:165];
        x.XTickLabel = [-500:50:320];
        x.XMinorTick = 'on';
        x.FontSize = txtsize;
        x.XTickLabelRotation = 90;
    
saveas(gcf, ['../Analysis Output/Poster_virtchan_GP.svg']);
close;
        
end
  

  
