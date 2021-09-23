%% Load conditions for subjects, average trials via ft_timelockanalysis

all_avg = struct();
gravg = struct();

%% skriv om och skriv bort sensspace och spara enligt formatet cond.XXlabel

for i = 1:length(sub_date.ID)
    
    subinpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    sensspace_avg.subjects{i,1} = sub_date.ID{i};
    
    for ii = 1%:length(conditions)
    
    tempdat = load([subinpath conditions{ii} 'ica.mat']);
    tempdat = tempdat.([conditions{ii} 'ica']);
    
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    trig = eval(['cond.' char(conditions(ii)) 'trig']);

    
        for iii = 5%nstim

        cfg = [];
        cfg.covariance = 'yes';
        cfg.covariancewindow = 'prestim';
        cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"?
        cfg.preproc.demean = 'yes';
        cfg.preproc.baselinewindow = [-0.200 -0.100]; %check this
        cfg.preproc.lpfilter = 'yes';
        cfg.preproc.lpfreq = 70;
        cfg.trials = tempdat.trialinfo == trig(iii);
        
        if ii == 1 & iii == 1
        PO60_70{i} = ft_timelockanalysis(cfg, tempdat);
        elseif ii == 1 & iii == 2
        PO60_75{i} = ft_timelockanalysis(cfg, tempdat);
        
        
        elseif ii == 1 & iii == 5
        PO60_90{i} = ft_timelockanalysis(cfg, tempdat);
        end
        
        %Collect all conditions averaged over trial for all subjects
        %sensspace_avg.(conditions{ii}){i, iii} = ft_timelockanalysis(cfg, tempdat);
        
        %Grab only 'avg' parameter in struct compatible with cat(mean())
        %all_avg.(conditions{ii}){i, iii} = sensspace_avg.(conditions{ii}){i, iii}.avg;
        
        %Calculate grand average over trial and subject per condition
        %gravg.(conditions{ii}){iii} = mean(cat(3, all_avg.(conditions{ii}){:, iii}), 3);
        
        %For trigs
        end
    
    %For conditions
    end

%For subjects
end

%save('../mat_data/sensspace/sensspace_avg.mat', 'sensspace_avg', '-v7.3');
%save('../mat_data/sensspace/all_avg.mat', 'all_avg', '-v7.3');
%save('../mat_data/sensspace/grand_avg.mat', 'gravg', '-v7.3');

%% Combine planar gradiometers

%% Butterfly with selected sensors

%Put grand average in a subjects structure for plotting
mean_sub = PO60_90{1,1};
mean_sub.avg = gravg.PO60{5};

cfg = [];
cfg.viewmode = 'butterfly';
cfg.channel = 'megmag';
cfg.colorgroups = 'allblack';
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.35 0];

ft_databrowser(cfg, mean_sub);

%% channels that we'll plot throughout

%Mesh from template brain
%Segment MRI
cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};

mri_segmented = ft_volumesegment(cfg, template_mri);

%Transform to neuromag coordsys - same as for sensors
mri_segmented = ft_convert_coordsys(mri_segmented, 'neuromag');

%Create mesh skull
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'skull';
cfg.numvertices = 2000;

mesh_skull = ft_prepare_mesh(cfg, mri_segmented);

ft_plot_mesh(mesh_skull);

%Plot higlighted sensors
colors = ones(306, 3);
l_chips = {'MEG1612', 'MEG1622', 'MEG1812', 'MEG1642', 'MEG1632', 'MEG1842', 'MEG1732', 'MEG1942', 'MEG1912'};
%r_chips = {'MEG2422', 'MEG2412', 'MEG2222', 'MEG2432', 'MEG2442', 'MEG2232', 'MEG2512', 'MEG2322', 'MEG2312'}

l_chips_idx = ismember(mean_sub.label, l_chips);

%Write colour vector to row in mean_sub.label (+1 to row for both grads on chip)
for i = 1:numel(mean_sub.label)
    
    if l_chips_idx(i) == 1
    colors(i, :) = [0 0 1];
    colors(i+1, :) = [0 0 1];
    end
    
end

figure('units', 'normalized', 'outerposition', [0 0 0.5 1]);
hold on
sensors = mean_sub.grad;
ft_plot_sens(sensors, 'facecolor', colors, 'facealpha', 0.8);

ft_plot_mesh(ft_convert_units(mesh_skull, 'cm'), 'edgecolor', [0.5 0.5 0.5]);

view([-100 25])

%MAGS
% cfg = [];
% cfg.parameter = 'avg';
% cfg.layout = 'neuromag306mag.lay';
% cfg.channel = 'megmag';
% cfg.showlabels = 'yes';
% 
% ft_multiplotER(cfg, mean_sub);

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
