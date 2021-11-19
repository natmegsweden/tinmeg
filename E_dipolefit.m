
%% Set up
%https://github.com/natmegsweden/meeg_course/blob/master/tutorial_04a_dipole_fitting.md

all_dip_pos = struct();

%Careful - Also used for dips2
toiN1 = [0.050 0.150];

%% Create initial individual dipole fit for PO and GO trials

for i = 1:length(sub_date.ID) %For all subjects
    
    %load headmodel and convert units
    headmodel = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/meg_headmodel.mat']);
    headmodel = headmodel.headmodel_meg;
    headmodel = ft_convert_units(headmodel, 'cm');
    
    %load MRI and convert units
    mri = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/mri_resliced.mat']);
    mri = mri.mri_resliced;
    mri = ft_convert_units(mri, 'cm');
    
    %Conditions PO60 and GO60
    for ii = [1, 5]
    
    %Set the stim trigger for PO or GO
    if ii == 1
        stim = 5
    elseif ii == 5
        stim = 1
    end

        for iii = stim
        
        %Load data
        ICAdat = load(['../mat_data/ICA/ID' sub_date.ID{i} '/' conditions{ii} 'ica.mat']);
        ICAdat = ICAdat.([conditions{ii} 'ica']);
        
        %Keep only MEG channels
        cfg = [];
        cfg.channel = 'MEG';
        ICAdat = ft_selectdata(cfg, ICAdat);
        
        % Whiten
        cfg = [];
        cfg.covariance          = 'yes';
        cfg.covariancewindow    = 'prestim';
        cfg.channel             = 'MEG';
        data_cov = ft_timelockanalysis(cfg, ICAdat);

        [u,s,v] = svd(data_cov.cov);
        d       = -diff(log10(diag(s)));
        d       = d./std(d);
        kappa   = find(d>5,1,'first');

        cfg            = [];
        cfg.channel    = 'meg';
        cfg.kappa      = kappa;
        dataw_meg      = ft_denoise_prewhiten(cfg, ICAdat, data_cov);
        
        clear ICAdat data_cov


        cfg = [];
        cfg.preproc.demean          = 'yes';
        cfg.preproc.baselinewindow  = [-0.200 0]; %150ms duration in sensspace processing?
        cfg.covariance              = 'yes';
        cfg.covariancewindow        = 'prestim';
        
        %cfg.preproc.lpfilter ???
        
        cfg.trials = dataw_meg.trialinfo == cond.([conditions{ii} 'trig'])(stim);
        
        timelockeds_wht = ft_timelockanalysis(cfg, dataw_meg);
        
        save(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/' cond.([conditions{ii} 'label']){stim} '_dip_timelockeds.mat'], 'timelockeds_wht')
        
        clear dataw_meg
        
        cfg = [];
        cfg.gridsearch          = 'yes';
        cfg.dipfit.metric       = 'rv';
        cfg.model               = 'regional';
        cfg.nonlinear           = 'yes';
        cfg.latency             = toiN1; % loads on top, needed for dips2
        cfg.symmetry            = [];
        cfg.headmodel           = headmodel;
        cfg.numdipoles          = 2;
        cfg.symmetry            = 'x';
        cfg.channel             = 'meg';
        cfg.resolution = 1;

        dips = ft_dipolefitting(cfg, timelockeds_wht);
        
        %Name for output files
        condlab = cond.([conditions{ii} 'label']){stim};
        
        %Plot subject dipole pos on mr and save for inspection
        figure; hold on
        pos = dips.dip.pos(1,:);

        ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', pos, 'orientation', [0 1 0]); hold on
        ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', pos, 'orientation', [1 0 0]); hold on
        ft_plot_dipole(dips.dip.pos(1,:), mean(dips.dip.mom(1:3,:),2), 'diameter', .5, 'unit', 'cm', 'color','g'); hold on
        ft_plot_dipole(dips.dip.pos(2,:), mean(dips.dip.mom(4:6,:),2), 'diameter', .5, 'unit', 'cm', 'color','g'); hold on
        view([0 1 0])
        
        saveas(gcf, ['../Analysis Output/dips/ID' sub_date.ID{i} '_' condlab '_dipoles.png']);
        close;
        
        %save dip.pos
        all_dip_pos.([condlab 'L'])(i,1:3) = dips.dip.pos(1,1:3);
        all_dip_pos.([condlab 'R'])(i,1:3) = dips.dip.pos(2,1:3);
        
        save(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/' condlab '_dipolefit.mat'], 'dips')
        
        clear dips timelockeds_wht
        
        end
        
    end
    
    clear d pos s u v kappa condlab headmodel mri stim i ii iii
    
end

%save(['../mat_data/source_reconstruction/all_dip_pos.mat'], 'all_dip_pos')

%% Regressing blinks/other processing errors

all_dip_pos = load('../mat_data/source_reconstruction/all_dip_pos.mat');
all_dip_pos = all_dip_pos.all_dip_pos;

%Subjects (idx) identified to have non-valid dipole positions (i.e. in eyes) - PO only!
for i = [4, 6, 16, 22];

    %Load dipolefit
    dips = load(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/PO60_90_dipolefit.mat']);
    dips = dips.dips;

    %load subject timelockeds, headmodel and mri
    headmodel = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/meg_headmodel.mat']);
    headmodel = headmodel.headmodel_meg;
    headmodel = ft_convert_units(headmodel, 'cm');

    mri = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/mri_resliced.mat']);
    mri = mri.mri_resliced;
    mri = ft_convert_units(mri, 'cm');

    timelockeds_wht = load(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/PO60_90_dip_timelockeds.mat']);
    timelockeds_wht = timelockeds_wht.timelockeds_wht;
    timelockeds_wht.time = round(timelockeds_wht.time, 4); %Round to use find()


    %Remove Vmodel from Vdata for subjects with dipoles fitted to blinks
    dif = dips.Vdata - dips.Vmodel;

    cropin = find(timelockeds_wht.time == toiN1(1));  %Find TOI start
    cropout = find(timelockeds_wht.time == toiN1(2)); %Find TOI end

    evoked_2 = timelockeds_wht; %Copy timelockeds
    evoked_2.avg(:, cropin:cropout) = dif; %Overwrite TOI with Vdata-Vmodel data

    %New dipolefit
    cfg = [];
    cfg.gridsearch          = 'yes';
    cfg.dipfit.metric       = 'rv';
    cfg.model               = 'regional';
    cfg.nonlinear           = 'yes';
    cfg.latency             = toiN1; %loads up top to be same as first dipfit
    cfg.symmetry            = [];
    cfg.headmodel           = headmodel;
    cfg.numdipoles          = 2;
    cfg.symmetry            = 'x';
    cfg.channel             = 'meg';
    cfg.resolution = 1;

    dips = ft_dipolefitting(cfg, evoked_2);
    
    %Plot bilat
    figure; hold on
    pos = dips.dip.pos(1,:);

    ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', pos, 'orientation', [0 1 0]); hold on
    ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', pos, 'orientation', [1 0 0]); hold on
    ft_plot_dipole(dips.dip.pos(1,:), mean(dips.dip.mom(1:3,:),2), 'diameter', .5, 'unit', 'cm', 'color','g'); hold on
    ft_plot_dipole(dips.dip.pos(2,:), mean(dips.dip.mom(4:6,:),2), 'diameter', .5, 'unit', 'cm', 'color','g'); hold on
    view([0 1 0])

    saveas(gcf, ['../Analysis Output/dips/ID' sub_date.ID{i} '_PO60_90_dipoles2.png']);
    close;

    %save dip.pos
    all_dip_pos.PO60_90L(i,1:3) = dips.dip.pos(1,1:3); %Static var name!
    all_dip_pos.PO60_90R(i,1:3) = dips.dip.pos(2,1:3); %Static var name!
    
    save(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/PO60_90_dipolefit2.mat'], 'dips')
    
end

clear timelockeds_wht pos mri i headmodel evoked_2 dips dif cropin cropout

%save(['../mat_data/source_reconstruction/all_dip_pos2.mat'], 'all_dip_pos')

%% Correct dip.pos to sourcemodel and create virtual channels

%To do: Split loop for correct dip.pos and ft_virtualchannel
%Include ft_timelockedanalysis before virtualchannel instead of loading whitened data from dips processing.
%Consider filtering as for sensspace processing

%NB - Load dip positions - 2 in filename is updated
all_dip_pos = load('../mat_data/source_reconstruction/all_dip_pos2.mat');
all_dip_pos = all_dip_pos.all_dip_pos;

%Empty struct for virtual channel data
all_virtch = struct();

for i = 1:length(sub_date.ID)

%load headmodel and convert units
headmodel = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/meg_headmodel.mat']);
headmodel = headmodel.headmodel_meg;
headmodel = ft_convert_units(headmodel, 'cm');

%load sourcemodel (created in E_beamform script)
source_org = load(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/source_org.mat']);
source_org = source_org.source_org;
source_org = ft_convert_units(source_org, 'cm');

%load MRI and convert units
mri = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/mri_resliced.mat']);
mri = mri.mri_resliced;
mri = ft_convert_units(mri, 'cm');
    
    for ii = [1, 5]%:length(conditions)
    
    nstim = length(eval(['cond.' conditions{ii} 'trig']));
    trig = eval(['cond.' conditions{ii} 'trig']);
    
    %Set the stim trigger for PO or GO
    if ii == 1
        stim = 5
    elseif ii == 5
        stim = 1
    end
    
    %Or load ICA and run ft_timlockedanalysis?
    %Loads whitened timelockeds frop dips
    timelockeds = load(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/' cond.([conditions{ii} 'label']){stim} '_dip_timelockeds.mat']);
    timelockeds = timelockeds.timelockeds_wht;

        for iii = stim
            
        %Warning re-using condlab
        condlab = cond.([conditions{ii} 'label']){stim};
        
            for side = ['L' 'R']; %Not really L/R but 1/2 is taken..

                %calculate corresponding closest source from dip.pos
                dist = nan(length(source_org.pos), 1);
                for jj = 1:length(source_org.pos)
                    dist(jj) =  norm(all_dip_pos.([condlab side])(i,:) - source_org.pos(jj,:));
                end

                [val, idx] = min(dist);

                all_dip_pos.([condlab side 'src'])(i,:) = source_org.pos(idx,:); %Save cordinates with dip.pos

                %Create virtual channel at dipole position
                cfg = [];
                cfg.pos = all_dip_pos.([condlab side 'src'])(i,:);
                cfg.method = 'svd';

                virtch = ft_virtualchannel(cfg, timelockeds, source_org)
                all_virtch.([condlab side])(i,:) = virtch.avg
                clear virtch

            %side (L/R);
            end
        
        %Trial
        end
       
    %Condition    
    end
    
    clear headmodel source_org mri
%Subject
end

clear trig val side nstim i ii iii jj idx condlab stim dist

save(['../mat_data/source_reconstruction/all_dip_pos3.mat'], 'all_dip_pos')
save(['../mat_data/source_reconstruction/all_virtch.mat'], 'all_virtch')

%% Plot dip.pos on template brain

%load data
load(['../mat_data/source_reconstruction/all_dip_pos.mat']);
load(['../mat_data/source_reconstruction/all_virtch.mat']);

load standard_mri

mri = ft_convert_units(mri, 'cm');

%Plot
figure; hold on
pos = all_dip_pos.PO60_90L(1,:);

ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', [0 0 0.5], 'orientation', [0 0 1]); hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', [0 0 0], 'orientation', [1 0 0]); hold on

for j = 1:22
ft_plot_dipole(all_dip_pos.PO60_90L(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','b', 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.PO60_90Lsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on


ft_plot_dipole(all_dip_pos.PO60_90R(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','r', 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.PO60_90Rsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on
end
view([0 0 1])

%% standard mesh brain with dip.pos and corrected (source) dip.pos for PO and GO

%Window/patch colors
red = [241 88 84]/256;
orange = [250 164 58]/256;
green = [96 189 104]/256;
blue = [93 165 218]/256;
purple = [178 118 178]/256;

%Standard headmodel from running standard_mri throught MRI_prep2 pipeline
load('../mat_data/MRI_mat/standard_headmodel.mat');

headmodel_std = ft_convert_units(headmodel_std, 'cm');

figure; subplot(1,2,1); hold on;
ft_plot_headmodel(headmodel_std, 'facealpha', 0.2, 'edgecolor', [0.5 0.5 0.5]);

%Puls
for j = 1:22
ft_plot_dipole(all_dip_pos.PO60_90L(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', blue, 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.PO60_90Lsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', red, 'alpha', 0.5); hold on

%for connecting lines "L"
x = [all_dip_pos.PO60_90Lsource(j,1) all_dip_pos.PO60_90L(j,1)];
y = [all_dip_pos.PO60_90Lsource(j,2) all_dip_pos.PO60_90L(j,2)];
z = [all_dip_pos.PO60_90Lsource(j,3) all_dip_pos.PO60_90L(j,3)];

plot3(x, y, z, 'Color', [0 0 0]);

%for connecting lines "R"
x = [all_dip_pos.PO60_90Rsource(j,1) all_dip_pos.PO60_90R(j,1)];
y = [all_dip_pos.PO60_90Rsource(j,2) all_dip_pos.PO60_90R(j,2)];
z = [all_dip_pos.PO60_90Rsource(j,3) all_dip_pos.PO60_90R(j,3)];

plot3(x, y, z, 'Color', [0 0 0]);

ft_plot_dipole(all_dip_pos.PO60_90R(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', blue, 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.PO60_90Rsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', red, 'alpha', 0.5); hold on
end
view([0.5 0.75 0])

%GAP
subplot(1,2,2); hold on;
ft_plot_headmodel(headmodel_std, 'facealpha', 0.2, 'edgecolor', [0.5 0.5 0.5]);

for j = 1:22
ft_plot_dipole(all_dip_pos.GO_60L(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', orange, 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.GO_60Lsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', red, 'alpha', 0.5); hold on

%for connecting lines "L"
x = [all_dip_pos.GO_60Lsource(j,1) all_dip_pos.GO_60L(j,1)];
y = [all_dip_pos.GO_60Lsource(j,2) all_dip_pos.GO_60L(j,2)];
z = [all_dip_pos.GO_60Lsource(j,3) all_dip_pos.GO_60L(j,3)];

plot3(x, y, z, 'Color', [0 0 0]);

%for connecting lines "R"
x = [all_dip_pos.GO_60Rsource(j,1) all_dip_pos.GO_60R(j,1)];
y = [all_dip_pos.GO_60Rsource(j,2) all_dip_pos.GO_60R(j,2)];
z = [all_dip_pos.GO_60Rsource(j,3) all_dip_pos.GO_60R(j,3)];

plot3(x, y, z, 'Color', [0 0 0]);

ft_plot_dipole(all_dip_pos.GO_60R(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', orange, 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.GO_60Rsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', red, 'alpha', 0.5); hold on
end
view([-0.5 0.75 0])


%% standard mesh brain with source dip.pos for PO and GO

figure; hold on;
ft_plot_headmodel(headmodel_std, 'facealpha', 0.2, 'edgecolor', [0.5 0.5 0.5]);

for j = 1:22
ft_plot_dipole(all_dip_pos.PO60_90Lsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', blue, 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.PO60_90Rsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', blue, 'alpha', 0.5); hold on  

ft_plot_dipole(all_dip_pos.GO_60Lsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', orange, 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.GO_60Rsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color', orange, 'alpha', 0.5); hold on
end
view([0.5 0.75 0])


%%

%Window/patch colors
red = [241 88 84]/256;
orange = [250 164 58]/256;
green = [96 189 104]/256;
blue = [93 165 218]/256;
purple = [178 118 178]/256;


% PULSE ONLY
figure('Position', [400 400 1800 400]); hold on;
xlim([41 165]);
%ylim([minylim maxylim]);
title('Pulse', 'Interpreter', 'none');
for i = 1:21%numel(sub_date.ID)
    
    plot(all_virtch.PO60_90R(i,:), 'Color', [0 0 0 0.5])
    tempmeanR(i,1:165) = all_virtch.PO60_90R(i,:);
    
    plot(all_virtch.PO60_90L(i,:), 'Color', [0 0 0 0.5])
    tempmeanL(i,1:165) = all_virtch.PO60_90L(i,:);
    
end

%N1 patch
patch('Faces', [1 2 3 4], 'Vertices', [111 -2; 111 2; 131 2; 131 -2], 'FaceColor', orange, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%P1 patch
%patch('Faces', [1 2 3 4], 'Vertices', [P2on minylim; P2on maxylim; P2off maxylim; P2off minylim], 'FaceColor', green, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

plot(mean(tempmeanR), 'Color', [0.75 0 0], 'LineWidth', 1.5); clear tempmean;
plot(mean(tempmeanL), 'Color', [0 0 0.75], 'LineWidth', 1.5); clear tempmean;
%plot([101 101], [minylim maxylim], 'k --');
x = gca;
x.XTick = [1:10:165];
x.XTickLabel = [-500:50:320];
x.XMinorTick = 'on';
x.FontSize = 20;
x.XTickLabelRotation = 90;


%GAP ONLY
figure('Position', [400 400 1800 400]); hold on;
xlim([41 165]);
%ylim([minylim maxylim]);
title('Gap', 'Interpreter', 'none');
for i = 1:21%numel(sub_date.ID)
    
    plot(all_virtch.GO_60R(i,:), 'Color', [0 0 0 0.5])
    tempmeanR(i,1:165) = all_virtch.GO_60R(i,:);
    
    plot(all_virtch.GO_60L(i,:), 'Color', [0 0 0 0.5])
    tempmeanL(i,1:165) = all_virtch.GO_60L(i,:);
    
end

%N1 patch
patch('Faces', [1 2 3 4], 'Vertices', [111 -2; 111 2; 131 2; 131 -2], 'FaceColor', orange, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

%P1 patch
%patch('Faces', [1 2 3 4], 'Vertices', [P2on minylim; P2on maxylim; P2off maxylim; P2off minylim], 'FaceColor', green, 'FaceAlpha', 0.1, 'EdgeAlpha', 0);

plot(mean(tempmeanR), 'Color', [0.75 0 0], 'LineWidth', 1.5); clear tempmean;
plot(mean(tempmeanL), 'Color', [0 0 0.75], 'LineWidth', 1.5); clear tempmean;
%plot([101 101], [minylim maxylim], 'k --');
x = gca;
x.XTick = [1:10:165];
x.XTickLabel = [-500:50:320];
x.XMinorTick = 'on';
x.FontSize = 20;
x.XTickLabelRotation = 90;
