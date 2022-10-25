%% To do:

% n = 8 PO amplitude response
% rank by inhibition

%% MultiplotER of inhib-vectors

bkgs = [0 3 8];

% For subjects
for i = 1:numel(sub_date.ID)

    for ii = 1:numel(bkgs)

    %load timelockeds
    PO = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/PO_0' num2str(bkgs(ii)) '_tlks_cmb.mat']);
    PO = PO.timelockeds_cmb;
    
    GPP = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/GPP_0' num2str(bkgs(ii)) '_tlks_cmb.mat']);
    GPP = GPP.timelockeds_cmb;
    
    %Select GRADS only
    cfg = [];
    cfg.channel = 'MEGGRAD'
    
    PO = ft_selectdata(cfg, PO);
    GPP = ft_selectdata(cfg, GPP);
    
    %Calculate inhib vector (PO-GPP)
    tempdat = PO.avg-GPP.avg;
    
    %Copy PO struct and overwrite with inhib vector to keep sensor/time/label fields
    tempstruct = PO;
    tempstruct.avg = tempdat;

    %Sort on max peak following pulse onset
    [val, ind] = sort(max(tempstruct.avg(:,101:end)'), 'descend');

    %Max MEAN amplitude following pulse
    %[val, ind] = sort(mean(tempstruct.avg(:,101:end), 2), 'descend') 
    
    %Colorgroups with ind2 as own group
    colgroups = repmat(1, 1, 102);
    colgroups(ind(1:10)) = 2;
    
    %MultiplotER
    cfg = [];
    cfg.parameter = 'avg';
    cfg.xlim = [-0.25 0.5];
    cfg.ylim = [-1 1];
    cfg.layout = 'neuromag306cmb.lay';
    cfg.colorgroups = colgroups;
    cfg.linecolor = [0 0 0; 1 0 0];
    ft_multiplotER(cfg, tempstruct)
    ax = gca;
    ax.Parent.Position = [500 100 1300 1100];

    title([sub_date.ID{i} ': Background: ' num2str(bkgs(ii))])

    %saveas(gcf, ['../Analysis Output/' [sub_date.ID{i} '_' num2str(bkgs(ii)) 'mplot.svg']]);
    %close

    %Plot highligted sensors together

    figure('Position', [400 400 1100 400]); hold on
    plot(tempstruct.avg(ind(1:10)',:)', 'Color', [0 0 0 0.9]);

    xlim([51 201])
    set(gca, 'YGrid', 'on')
    ylim([-1 1]);
    xline([101])
    xticks(1:25:201)
    xticklabels(-500:125:500)
    yline([0], 'Color', [1 0 0])

    title([sub_date.ID{i} ': Background: ' num2str(bkgs(ii)) ', soi']);

    saveas(gcf, ['../Analysis Output/' [sub_date.ID{i} '_' num2str(bkgs(ii)) '_soi.svg']]);
    close

    clear tempdat

    %For bkgs
    end

    %For subjects
end

%% Inhib-vector plots

bkgs = [0 3 8];

% For subjects
for i = 1:numel(sub_date.ID)

    for ii = 1:numel(bkgs)

    %load timelockeds
    PO = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/PO_0' num2str(bkgs(ii)) '_tlks_cmb.mat']);
    PO = PO.timelockeds_cmb;
    
    GPP = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/GPP_0' num2str(bkgs(ii)) '_tlks_cmb.mat']);
    GPP = GPP.timelockeds_cmb;
    
    %Select GRADS only
    cfg = [];
    cfg.channel = 'MEGGRAD'
    
    PO = ft_selectdata(cfg, PO);
    GPP = ft_selectdata(cfg, GPP);
    
    %Calculate inhib vector (PO-GPP)
    tempdat = PO.avg-GPP.avg;
    
    %Plot parameters
    patch_range = 4;
    %rois = [119 135 150];
    ylims = [-8*10^-12 8*10^-12];

    %Plot inihibition amplitude for all sensors    
    figure('Position', [400 400 1100 400]); hold on
    plot(tempdat', 'Color', [0 0 0 0.2])
    plot(mean(tempdat', 2), 'Color', [1 0 0 0.75], 'LineWidth', 1.5)
    xlim([51 201])
    set(gca, 'YGrid', 'on')
    ylim(ylims)
    %xlim([25 150])
    xline([101])
    xticks(1:25:201)
    xticklabels(-500:125:500)
    yline([0], 'Color', [1 0 0])

    title([sub_date.ID{i} ': Background: ' num2str(bkgs(ii))])

    saveas(gcf, ['../Analysis Output/' [sub_date.ID{i} '_' num2str(bkgs(ii)) '.svg']]);
    close
    clear tempdat
    
    %patch([rois(1)-patch_range rois(1)-patch_range rois(1)+patch_range rois(1)+patch_range], [ylims(1) ylims(2) ylims(2) ylims(1) ], [0 0 1], 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
    %patch([rois(2)-patch_range rois(2)-patch_range rois(2)+patch_range rois(2)+patch_range], [ylims(1) ylims(2) ylims(2) ylims(1) ], [0 0 1], 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
    %patch([rois(3)-patch_range rois(3)-patch_range rois(3)+patch_range rois(3)+patch_range], [ylims(1) ylims(2) ylims(2) ylims(1) ], [0 0 1], 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
    
    %For bkgs ii
    end

%For subjects i
end

%Identify different components (in time)

%Heatmap on helmet

%Highest L/R