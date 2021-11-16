
%https://github.com/natmegsweden/meeg_course/blob/master/tutorial_04a_dipole_fitting.md

all_dip_pos = struct();
all_virtch = struct();

for i = 1:length(sub_date.ID)

    %load leadfield (from E_beamformer) - Subject to change
    %leadfield = load(['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/leadfield.mat']);
    %leadfield = leadfield.leadfield;
    
    %load headmodel (from D_MR_prep2)
    headmodel = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/meg_headmodel.mat']);
    headmodel = headmodel.headmodel_meg;
    
    % Convert units
    headmodel = ft_convert_units(headmodel, 'cm');
    
    %ft_plot_headmodel(headmodel);
    
    %Load sourcemodel (created in E_beamform script)
    load(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/source_org.mat']);
    source_org = ft_convert_units(source_org, 'cm')
    
    %load reslice MRI for plot
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

        for iii = stim
        
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
       
        %implement adaptive baseline window and trigger
   
        cfg = [];
        cfg.preproc.demean          = 'yes';
        cfg.preproc.baselinewindow  = [-0.200 0];
        cfg.covariance              = 'yes';
        cfg.covariancewindow        = 'prestim';
        
        cfg.trials = dataw_meg.trialinfo == cond.([conditions{ii} 'trig'])(stim);
        
        evoked_wht = ft_timelockanalysis(cfg, dataw_meg);
        
        save(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/' cond.([conditions{ii} 'label']){stim} 'dip_timelockeds.mat'], 'evoked_wht')
        
        clear dataw_meg
        
        toiN1 = [0.050 0.150];
        
        cfg = [];
        cfg.gridsearch          = 'yes';
        cfg.dipfit.metric       = 'rv';
        cfg.model               = 'regional';
        cfg.nonlinear           = 'yes';
        cfg.latency             = toiN1; % FIRST COMPONENT
        cfg.symmetry            = [];
        cfg.headmodel           = headmodel;
        cfg.numdipoles          = 2;              % we expect bilateral activity
        cfg.symmetry            = 'x';
        cfg.channel             = 'meg';
        cfg.resolution = 1;
                
        dips = ft_dipolefitting(cfg, evoked_wht);
                
        %Convert units for plotting
        %dip_all_bil_wht.dip = ft_convert_units(dip_all_bil_wht.dip, 'mm');

        %Create name struct to maintain some readability
        name{1} = ([cond.([conditions{ii} 'label']){stim} 'L']);
        name{2} = ([cond.([conditions{ii} 'label']){stim} 'R']);
        name{3} = cond.([conditions{ii} 'label']){stim};
        name{4} = ([cond.([conditions{ii} 'label']){stim} 'Lsource']);
        name{5} = ([cond.([conditions{ii} 'label']){stim} 'Rsource']);
        
        %Plot bilat
        figure; hold on
        pos = dips.dip.pos(1,:);

        ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', pos, 'orientation', [0 1 0]); hold on
        ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', pos, 'orientation', [1 0 0]); hold on
        ft_plot_dipole(dips.dip.pos(1,:), mean(dips.dip.mom(1:3,:),2), 'diameter', .5, 'unit', 'cm', 'color','g'); hold on
        ft_plot_dipole(dips.dip.pos(2,:), mean(dips.dip.mom(4:6,:),2), 'diameter', .5, 'unit', 'cm', 'color','g'); hold on
        view([0 1 0])
        
        saveas(gcf, ['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/' name{3} '_dipoles.png']);
        close;

        
        %save dip.pos
        all_dip_pos.(name{1})(i,1:3) = dips.dip.pos(1,1:3);
        all_dip_pos.(name{2})(i,1:3) = dips.dip.pos(2,1:3);
        
        save(['../mat_data/source_reconstruction/ID' sub_date.ID{i} '/' name{3} '_dipolefit.mat'], 'dips')
        
        clear dips
        
        %calculate corresponding closest source from dip.pos
        %Left
        dist = nan(length(source_org.pos), 1);
        for jj = 1:length(source_org.pos)
            dist(jj) =  norm(all_dip_pos.(name{1})(i,:) - source_org.pos(jj,:));
        end

        [val, idx] = min(dist);

        all_dip_pos.(name{4})(i,:) = source_org.pos(idx,:); %Save cordinates with dip.pos
        
        %Right
        dist = nan(length(source_org.pos), 1);
        for jj = 1:length(source_org.pos)
            dist(jj) =  norm(all_dip_pos.(name{2})(i,:) - source_org.pos(jj,:));
        end

        [val, idx] = min(dist);

        all_dip_pos.(name{5})(i,:) = source_org.pos(idx,:); %Save cordinates with dip.pos
        
        %Create virtual channel at dipole position
        cfg = [];
        cfg.pos = all_dip_pos.(name{5})(i,:);
        cfg.method = 'svd';
        
        virtchR = ft_virtualchannel(cfg, evoked_wht, source_org)
        all_virtch.([name{2}])(i,:) = virtchR.avg
        clear virtchR
        
        cfg.pos = all_dip_pos.(name{4})(i,:);
        virtchL = ft_virtualchannel(cfg, evoked_wht, source_org)
        all_virtch.([name{1}])(i,:) = virtchL.avg
        clear virtchL
        
        end
        
    end
    
end

save(['../mat_data/source_reconstruction/all_dip_pos.mat'], 'all_dip_pos')
save(['../mat_data/source_reconstruction/all_virtch.mat'], 'all_virtch')

%% Regressing blinks/other processing errors

%To do: load identified subjects, rename dips.mat for book-keeping and write new dips.mat below

% Blinks:
% 0756 PO - 4 - dips2 ok
% 0835 PO - 6 - dips2 ok
% 0863 PO - 22 - dips2 ok

% Other:
% 0853 PO outside above/motor? - 16 - dips2 ok

dips = load(['../mat_data/source_reconstruction/ID0853/PO60_90_dipolefit.mat']);
dips = dips.dips;

%load subject evoked_wht and headmodel

load standard_mri;
mri = ft_convert_units(mri, 'cm');

%Plot
figure; hold on
pos = dips.dip.pos(1,:);
ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', [0 0 0.5], 'orientation', [0 0 1]); hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', [0 0 0], 'orientation', [1 0 0]); hold on

ft_plot_dipole(dips.dip.pos(1,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on
ft_plot_dipole(dips.dip.pos(2,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on

%Remove Vmodel from Vdata for subjects with dipoles fitted to blinks
dif = dips.Vdata - dips.Vmodel;

evoked_2 = evoked_wht; %Copy timelockeds
evoked_2.avg(:, 111:131) = dif; %Overwrite TOI with Vdata-Vmodel data

%New dipolefit
toiN1 = [0.050 0.150];

cfg = [];
cfg.gridsearch          = 'yes';
cfg.dipfit.metric       = 'rv';
cfg.model               = 'regional';
cfg.nonlinear           = 'yes';
cfg.latency             = toiN1; % FIRST COMPONENT
cfg.symmetry            = [];
cfg.headmodel           = headmodel;
cfg.numdipoles          = 2;              % we expect bilateral activity
cfg.symmetry            = 'x';
cfg.channel             = 'meg';
cfg.resolution = 1;

dips2 = ft_dipolefitting(cfg, evoked_2);

%Plot 2
figure; hold on
pos = dips2.dip.pos(1,:);
ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', [0 0 0.5], 'orientation', [0 0 1]); hold on
ft_plot_slice(mri.anatomy, 'transform', mri.transform,'location', [0 0 0], 'orientation', [1 0 0]); hold on

ft_plot_dipole(dips2.dip.pos(1,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on
ft_plot_dipole(dips2.dip.pos(2,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on



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

%% Mesh brain

%Standard headmodel from running standard_mri throught MRI_prep2 pipeline
load('../mat_data/MRI_mat/standard_headmodel.mat');

headmodel_std = ft_convert_units(headmodel_std, 'cm');

figure; hold on;
ft_plot_headmodel(headmodel_std, 'facealpha', 0.5, 'edgecolor', [0.5 0.5 0.5]);

for j = 1:21
ft_plot_dipole(all_dip_pos.PO60_90L(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','b', 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.PO60_90Lsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on


ft_plot_dipole(all_dip_pos.PO60_90R(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','r', 'alpha', 0.5); hold on
ft_plot_dipole(all_dip_pos.PO60_90Rsource(j,:), [0 0 0], 'diameter', 0.5, 'unit', 'cm', 'color','g', 'alpha', 0.5); hold on
end




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
