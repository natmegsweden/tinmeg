%% Load vars and gathered timelocked data

% Run setup for common variables if not already loaded
if exist('sub_date', 'var') == 0; run Conditions_triggers.m; end;

tlk_sub_cmb = load('../mat_data/timelockeds/tlk_sub_cmb.mat');
tlk_sub_cmb = tlk_sub_cmb.tlk_sub_cmb;

%Load test data from one example subject
timelockeds_cmb = load(['../mat_data/timelockeds/ID' sub_date.ID{1} '/PO60_90_tlks_cmb.mat']);
timelockeds_cmb = timelockeds_cmb.timelockeds_cmb;

timelockeds = load(['../mat_data/timelockeds/ID' sub_date.ID{1} '/PO60_90_tlks.mat']);
timelockeds = timelockeds.timelockeds;

%Get sensor positions and labels
sensors = timelockeds.grad;

senslab = timelockeds_cmb.label;
senspos = timelockeds_cmb.grad.chanpos;
timevec = round(timelockeds_cmb.time,3);

clear timlockeds_cmb timelockeds

%To do
%TOIon TOIoff vars

%% Sensshape plot with L and R chips highlighted

if ~exist('mri_segmented', 'var')
    
    load standard_mri;
    template_mri = mri;
    clear mri;
    
    cfg = [];
    cfg.output = 'scalp';
    
    %Segment template MRI for head mesh
    mri_segmented = ft_volumesegment(cfg, template_mri);
    
    %Transform to neuromag coordsys - same as for sensors
    mri_segmented = ft_convert_coordsys(mri_segmented, 'neuromag');
end

if ~exist('mesh_scalp', 'var')
    %Create mesh skull
    cfg = [];
    cfg.method = 'projectmesh';
    cfg.tissue = 'scalp';
    cfg.numvertices = 1000;
    
    mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);
    
    %Move scalp to not clip through sensors and be positioned nicely in helmet
    mesh_scalp.pos(:,3) = mesh_scalp.pos(:,3) - 35; %Down
    mesh_scalp.pos(:,2) = mesh_scalp.pos(:,2) - 15; %Back
end

%Create color matrix
colors2 = zeros(306, 3);

colors2(:,1) = sensors.chanpos(:,1) > 0;
colors2(:,3) = sensors.chanpos(:,1) < 0;

figure('Position', [400 200 1800 1000]); hold on
subplot(1,2,1)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

subplot(1,2,2)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([-100 25])

%% Find highest POXX_90-response grad for Left and Right hemishphere per subject

%Prepare structures
L_topgrads = struct();
R_topgrads = struct();

%Define TOI 1 (50 - 150 ms)
toi1 = find(timevec == 0.050);
toi2 = find(timevec == 0.150);

%Get labels for Right and Left sensors
Rchan_lab = senslab(senspos(:,1) > 0, :);
Lchan_lab = senslab(senspos(:,1) < 0, :);

for ii = [1, 2] %For condition PO60 and PO70 (index in var: 'conditions')

    for i = 1:numel(sub_date.ID)
        
        %Select stim level 90dB as index in var: 'cond'
        if ii == 1
            iii = 5;
        elseif ii == 2;
            iii = 4;
        end

        %The condition and stimuli
        %name = cond.([conditions{ii} 'label']){iii};

        %Define Right and Left sensors based on x-coordinate of sensor position
        Rchan = tlk_sub_cmb.(conditions{ii}){i, iii}(senspos(:,1) > 0, :);
        Lchan = tlk_sub_cmb.(conditions{ii}){i, iii}(senspos(:,1) < 0, :);
        
        %Find sensor on RIGHT side with biggest amplitude response in TOI
        [datR, indR] = sort(mean(Rchan(:,toi1:toi2), 2), 'descend');
        topRname = Rchan_lab{indR(1)};
        topRind = find(ismember(senslab, topRname));
        
        %Find sensor on LEFT side with biggest amplitude response in TOI
        [datL, indL] = sort(mean(Lchan(:,toi1:toi2), 2), 'descend');
        topLname = Lchan_lab{indL(1)};
        topLind = find(ismember(senslab, topLname));
    
        %Write top channel index (out of 204)
        L_topgrads.([conditions{ii} 'chind']){i,:} = topLind;
        R_topgrads.([conditions{ii} 'chind']){i,:} = topRind;

        %Write top channel label/name
        L_topgrads.([conditions{ii} 'chan']){i,:} = topLname;
        R_topgrads.([conditions{ii} 'chan']){i,:} = topRname;

    %For subject
    end

%For conditions
end

clear topLname topRname topLind topRind;

%Count up what sensors are most common as top response
[count60L, name60L] = groupcounts(L_topgrads.PO60chan);
[count60R, name60R] = groupcounts(R_topgrads.PO60chan);

[count70L, name70L] = groupcounts(L_topgrads.PO70chan);
[count70R, name70R] = groupcounts(R_topgrads.PO70chan);

%Find and print info on most common top gradiometer
[C, I] = sort(count60L, 'descend');
top60L = name60L{I(1)};
disp(['Top sensor LEFT in 60 dB carrier is: ' top60L ' - highest response in ' num2str(count60L(I(1))) ' of ' num2str(numel(sub_date.ID)) ' participants']);

[C, I] = sort(count60R, 'descend');
top60R = name60R{I(1)};
disp(['Top sensor RIGHT in 60 dB carrier is: ' top60R ' - highest response in ' num2str(count60R(I(1))) ' of ' num2str(numel(sub_date.ID)) ' participants']);

[C, I] = sort(count70L, 'descend');
top70L = name70L{I(1)};
disp(['Top sensor LEFT in 70 dB carrier is: ' top70L ' - highest response in ' num2str(count70L(I(1))) ' of ' num2str(numel(sub_date.ID)) ' participants']);

[C, I] = sort(count70R, 'descend');
top70R = name70R{I(1)};
disp(['Top sensor RIGHT in 70 dB carrier is: ' top70R ' - highest response in ' num2str(count70R(I(1))) ' of ' num2str(numel(sub_date.ID)) ' participants']);

clear count60L name60L count70L name70R;

%% Gather response from all conditions from SOIs identified for POXX_90

topgrad_dat = struct();

for ii = 1:numel(conditions); %For condition PO60 and GP60 (index [1,3] in var: 'conditions')

    for iii = 1:numel(cond.([conditions{ii} 'label']));
    
            %The condition and stimuli
            name = cond.([conditions{ii} 'label']){iii};
            disp(name); %helps to monitor what's going on
    
        for i = 1:numel(sub_date.ID)
            
            if ii == 1 || ii == 3 || ii == 5 && iii == 1 %if 60 dB carrier
                
                  % Use most common top sensor
                  temp = tlk_sub_cmb.(conditions{ii}){i,iii};
                  topgrad_dat.([conditions{ii} '_L']){1,iii}(i,:) = temp(find(ismember(senslab, top60L)),:); 
                  topgrad_dat.([conditions{ii} '_R']){1,iii}(i,:) = temp(find(ismember(senslab, top60R)),:);
                  clear temp

%                 % This snippet is for unique top sensor per participant
%                 temp = tlk_sub_cmb.(conditions{ii}){i,iii};
%                 topgrad_dat.([conditions{ii} '_L']){1,iii}(i,:) = temp(L_topgrads.PO60chind{i},:); 
% 
%                 temp = tlk_sub_cmb.(conditions{ii}){i,iii};
%                 topgrad_dat.([conditions{ii} '_R']){1,iii}(i,:) = temp(R_topgrads.PO60chind{i},:);

            elseif ii == 2 || ii == 4 || ii == 5 && iii == 2 %if 70 dB carrier

                  % Use most common top sensor
                  temp = tlk_sub_cmb.(conditions{ii}){i,iii};
                  topgrad_dat.([conditions{ii} '_L']){1,iii}(i,:) = temp(find(ismember(senslab, top70L)),:); 
                  topgrad_dat.([conditions{ii} '_R']){1,iii}(i,:) = temp(find(ismember(senslab, top70R)),:);
                  clear temp
                
%                 % This snippet is for unique top sensor per participant
%                 temp = tlk_sub_cmb.(conditions{ii}){i,iii};
%                 topgrad_dat.([conditions{ii} '_L']){1,iii}(i,:) = temp(L_topgrads.PO70chind{i},:); 
% 
%                 temp = tlk_sub_cmb.(conditions{ii}){i,iii};
%                 topgrad_dat.([conditions{ii} '_R']){1,iii}(i,:) = temp(R_topgrads.PO70chind{i},:);

            end
        %For subject
        end

    %For stims 
    end

%For conditions
end


%% Figures of top gradiometers

%patch 112 132

% [185, 202, 254] light blue RGB values
% [1, 55, 203] dark blue RGB values
% [241, 215, 177] light orange RGB values
% [252, 96, 10] dark orange RGB values

%Blue color gradient
bluecol(:,1) = linspace(185, 1, 6) ./256; % R
bluecol(:,2) = linspace(202, 55, 6) ./256; % G
bluecol(:,3) = linspace(254, 203, 6) ./256; % B

%Orange color gradient
orgcol(:,1) = linspace(241, 252, 6) ./256; % R
orgcol(:,2) = linspace(215, 96, 6) ./256; % G
orgcol(:,3) = linspace(177, 10, 6) ./256; % B

%PO60
figure('Position', [774,254,745,832], 'Renderer','painters'); tiledlayout(3,1, 'TileSpacing','tight', 'TileIndexing','columnmajor');
nexttile; hold on;

patch([112 112 132 132], [0 1.3*10^-11 1.3*10^-11 0], [0.2 0.2 0.2], 'FaceAlpha', 0.125, 'EdgeAlpha', 0)

for i = 1:6 % Six stim (0-240ms ISI) in GP60 condition
    plot(mean(topgrad_dat.PO60_L{1,i}', 2), 'Color', bluecol(i,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.PO60_R{1,i}', 2), 'Color', orgcol(i,:), 'LineWidth', 1.5)
end

%legend({}, 'Location', 'northwest');

ylabel('ERF amplitude (T/cm)');

xline(101, '--k')

% xline(111, 'k')
% xline(131, 'k')

xline(101, '--k')
xticks([1:20:165])
xticklabels([])
xlim([41 161])

ylim([0 1.3*10^-11])

ax = gca;
ax.XGrid = 'on';

%saveas(gcf, '../Analysis Output/L_R_POtest.svg')


%GP60
% Consider line() to mark gap locations

yoffset = 0.25*10^-11;

nexttile; hold on;

patch([112 112 132 132], [0 1.3*10^-11 1.3*10^-11 0], [0.2 0.2 0.2], 'FaceAlpha', 0.125, 'EdgeAlpha', 0)

for i = 1:4; % Four stim (ISI0-240) in GP60 condition flip for order preference
    if i == 1;
    plot(mean(topgrad_dat.GP60_L{1,i}', 2), 'Color', bluecol(6,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.GP60_R{1,i}', 2), 'Color', orgcol(6,:), 'LineWidth', 1.5)
    else
    plot(mean(topgrad_dat.GP60_L{1,i}', 2)+(i-1)*yoffset, 'Color', bluecol(6,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.GP60_R{1,i}', 2)+(i-1)*yoffset, 'Color', orgcol(6,:), 'LineWidth', 1.5)
    end
end

ylim([0 1.3*10^-11])

xline(101, '--k')
xticks([1:20:165])
xticklabels([])
xlim([41 161])

ylabel('Relative ERF amplitude');

ax = gca;
ax.XGrid = 'on';

plot([43 53 67 77 79 89 91 101], [0.9 0.9 0.61 0.61 0.38 0.38 0.1 0.1]*10^-11, '|k', 'MarkerSize', 22)

nexttile; hold on;

patch([112 112 132 132], [0 1.3*10^-11 1.3*10^-11 0], [0.2 0.2 0.2], 'FaceAlpha', 0.125, 'EdgeAlpha', 0)

plot(mean(topgrad_dat.PO60_R{1,5}', 2), 'Color', orgcol(6,:), 'LineWidth', 1.5)
plot(mean(topgrad_dat.PO60_L{1,5}', 2), 'Color', bluecol(6,:), 'LineWidth', 1.5)

%GPISI240
plot(mean(topgrad_dat.GP60_R{1,4}', 2), ':', 'Color', orgcol(6,:), 'LineWidth', 1.5)
plot(mean(topgrad_dat.GP60_L{1,4}', 2), ':', 'Color', bluecol(6,:), 'LineWidth', 1.5)

legend({'', 'Pulse only RIGHT', 'Pulse only LEFT', 'Gap+Pulse (isi 240ms) RIGHT', 'Gap+Pulse (isi 240ms) LEFT'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

ylim([0 1.3*10^-11])

xline(101, '--k')
xticks([1:20:165])
xticklabels([-500:100:320])
xlim([41 161])

ylabel('ERF amplitude (T/cm)');
xlabel('Time (ms relative pulse onset)')

ax = gca;
ax.XGrid = 'on';

%saveas(gcf, '../Analysis Output/Figure1amps.svg');
%close;

%% Sensshape with gradients for n of top sensor

if ~exist('mri_segmented', 'var')
    
    load standard_mri;
    template_mri = mri;
    clear mri;
    
    cfg = [];
    cfg.output = 'scalp';
    
    %Segment template MRI for head mesh
    mri_segmented = ft_volumesegment(cfg, template_mri);
    
    %Transform to neuromag coordsys - same as for sensors
    mri_segmented = ft_convert_coordsys(mri_segmented, 'neuromag');
end

if ~exist('mesh_scalp', 'var')
    %Create mesh skull
    cfg = [];
    cfg.method = 'projectmesh';
    cfg.tissue = 'scalp';
    cfg.numvertices = 1000;
    
    mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);
    
    %Move scalp to not clip through sensors and be positioned nicely in helmet
    mesh_scalp.pos(:,3) = mesh_scalp.pos(:,3) - 35; %Down
    mesh_scalp.pos(:,2) = mesh_scalp.pos(:,2) - 15; %Back
end

%Create color matrix
colors2 = ones(306, 3);

%Blue color gradient
bluecol(:,1) = linspace(185, 1, 6) ./256; % R
bluecol(:,2) = linspace(202, 55, 6) ./256; % G
bluecol(:,3) = linspace(254, 203, 6) ./256; % B

%Orange color gradient
orgcol(:,1) = linspace(241, 252, 6) ./256; % R
orgcol(:,2) = linspace(215, 96, 6) ./256; % G
orgcol(:,3) = linspace(177, 10, 6) ./256; % B

%Color left blue
%Find sensor label and put in colorspace
sensind = find(ismember(sensors.label, top60L(1:7))); %Search for first gradiometer
colors2(sensind,:) = bluecol(6,:); %Color first gradiometer
colors2(sensind-1,:) = bluecol(6,:); %Color other magnetometer
colors2(sensind+1,:) = bluecol(6,:); %Color other gradiomter

%Color right orange
%Find sensor label and put in colorspace
sensind = find(ismember(sensors.label, top60R(1:7))); %Search for first gradiometer
colors2(sensind,:) = orgcol(6,:); %Color first gradiometer
colors2(sensind-1,:) = orgcol(6,:); %Color other magnetometer
colors2(sensind+1,:) = orgcol(6,:); %Color other gradiomter

figure('Position', [400 200 1800 1000], 'Renderer','painters'); hold on
subplot(1,2,1)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.9);
%ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

subplot(1,2,2)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.9);
%ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([-100 25])

saveas(gcf, '../Analysis Output/headshape_L_R.svg');
close;

%% TopoplotER

%load one subject data for approriate structure
meansub = load(['../mat_data/timelockeds/ID' num2str(sub_date.ID{1}) '/PO60_90_tlks_cmb.mat']);
meansub = meansub.timelockeds_cmb;

zlimlow = -1*10^-13;
zlimhigh = 1.5*10^-13;

%%%%%%%%%%%%%%%%%%%%%%
%PO60_90
%concatenate PO60_90 arrays
PO60_90_ga = cat(3,tlk_sub_cmb.PO60{:,5});
%grand average of all subjects
PO60_90_ga = mean(PO60_90_ga, 3);

%Put grand average array in meansub avg
meansub.avg = PO60_90_ga;

figure('Position', [200 200 800 800], 'Renderer','painters');
cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

cfg.xlim = [0.055 0.155];
cfg.baseline = [-0.150 0];

ft_topoplotER(cfg, meansub);

%title((['90dB pulse only 50-150ms']), 'Interpreter', 'none');
%colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
colormap(viridis) % change the colormap
colorbar

saveas(gcf, '../Analysis Output/topo_PO6090_N1.svg');
close;

%%%%%%%%%%%%%%%%%%%%%%
%GP60_i240
%concatenate PO60_90 arrays
GP60_i240_ga = cat(3,tlk_sub_cmb.GP60{:,4});
%grand average of all subjects
GP60_i240_ga = mean(GP60_i240_ga, 3);

%Put grand average array in meansub avg
meansub.avg = GP60_i240_ga;

figure('Position', [200 200 800 800], 'Renderer','painters');
cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'neuromag306mag.lay';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.zlim = [zlimlow zlimhigh];

cfg.xlim = [0.055 0.155];
cfg.baseline = [-0.440 -0.290];

ft_topoplotER(cfg, meansub);

%title((['Gap+Pulse 50-150ms (ISI 240ms)']), 'Interpreter', 'none');
%colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
colormap(viridis) % change the colormap
colorbar

saveas(gcf, '../Analysis Output/topo_GP60i240_N1.svg');
close;

%% Quantify response - Gather mean in TOI

%Gather mean amplitude in TOI
amp = struct();

topgrad_dat.PO60_L{1,1}(:, 112:132)
sides = ['L', 'R'];

%For side
for i = 1:numel(sides)

    %for conditions
    for ii = 1:numel(conditions)
    
        %for stim
        for iii = 1:numel(cond.([conditions{ii} 'label']))

            amp.([conditions{ii} '_' sides(i)])(:,iii) = mean(topgrad_dat.([conditions{ii} '_' sides(i)]){1,iii}(:, 112:132),2); %112:132 = 55-155ms

        %for stim
        end

    %for conditions
    end

%for side
end

csvwrite('../R data/amp_PO60_L', amp.PO60_L);
csvwrite('../R data/amp_PO70_L', amp.PO70_L);
csvwrite('../R data/amp_GP60_L', amp.GP60_L);
csvwrite('../R data/amp_GP70_L', amp.GP70_L);
csvwrite('../R data/amp_GO_L', amp.GO_L);

csvwrite('../R data/amp_PO60_R', amp.PO60_R);
csvwrite('../R data/amp_PO70_R', amp.PO70_R);
csvwrite('../R data/amp_GP60_R', amp.GP60_R);
csvwrite('../R data/amp_GP70_R', amp.GP70_R);
csvwrite('../R data/amp_GO_R', amp.GO_R);

%% Boxplots

blue = [0.0039 0.2148 0.7929];
orange = [0.9844 0.3750 0.0391];

figure('Renderer','painters'); tiledlayout(2,2, 'TileSpacing','tight');
nexttile; hold on;
boxchart(amp.PO60_L(:,:), 'BoxFaceColor', blue, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', blue, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'70', '75', '80', '85', '90', '95'})
ylabel('N1 ERF amplitude (T/cm)')
xlabel('Pulse level (dB)')

nexttile; hold on;
boxchart(amp.PO60_R(:,:), 'BoxFaceColor', orange, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', orange, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'70', '75', '80', '85', '90', '95'})
xlabel('Pulse level (dB)')


nexttile; hold on;
boxchart(amp.GP60_L(:,:), 'BoxFaceColor', blue, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', blue, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'0', '60', '120', '240'})
ylabel('N1 ERF amplitude (T/cm)')
xlabel('ISI (ms)')

nexttile; hold on;
boxchart(amp.GP60_R(:,:), 'BoxFaceColor', orange, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', orange, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'0', '60', '120', '240'})
xlabel('ISI (ms)')


%% Form 9 Fig 2

%patch 112 132

% [185, 202, 254] light blue RGB values
% [1, 55, 203] dark blue RGB values
% [241, 215, 177] light orange RGB values
% [252, 96, 10] dark orange RGB values

txtsize = 12;

%Blue color gradient
bluecol(:,1) = linspace(185, 1, 6) ./256; % R
bluecol(:,2) = linspace(202, 55, 6) ./256; % G
bluecol(:,3) = linspace(254, 203, 6) ./256; % B

%Orange color gradient
orgcol(:,1) = linspace(241, 252, 6) ./256; % R
orgcol(:,2) = linspace(215, 96, 6) ./256; % G
orgcol(:,3) = linspace(177, 10, 6) ./256; % B

%PO60
figure('Units', 'centimeters', 'Position',  [5 5 32 20], 'Renderer','painters');; tiledlayout(3,4, 'TileSpacing','tight', 'TileIndexing','columnmajor');
nexttile(1,1:2); hold on;

patch([112 112 132 132], [0 1.3*10^-11 1.3*10^-11 0], [0.2 0.2 0.2], 'FaceAlpha', 0.125, 'EdgeAlpha', 0)

for i = 1:6 % Six stim (0-240ms ISI) in GP60 condition
    plot(mean(topgrad_dat.PO60_L{1,i}', 2), 'Color', bluecol(i,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.PO60_R{1,i}', 2), 'Color', orgcol(i,:), 'LineWidth', 1.5)
end

%legend({}, 'Location', 'northwest');

set(gca, 'FontSize', txtsize)
ylabel('ERF amp (T/cm)');

xline(101, '--k')

% xline(111, 'k')
% xline(131, 'k')

xline(101, '--k')
xticks([1:20:165])
xticklabels([])
xlim([41 161])

ylim([0 1.3*10^-11])

ax = gca;
ax.XGrid = 'on';

%Fake legend to modify
legend({'', '70', '75', '80', '85', '90', '95', 'LEFT', 'RIGHT'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

%saveas(gcf, '../Analysis Output/L_R_POtest.svg')


%GP60
% Consider line() to mark gap locations

yoffset = 0.25*10^-11;

nexttile(2,1:2); hold on;

patch([112 112 132 132], [0 1.3*10^-11 1.3*10^-11 0], [0.2 0.2 0.2], 'FaceAlpha', 0.125, 'EdgeAlpha', 0)

for i = 1:4; % Four stim (ISI0-240) in GP60 condition flip for order preference
    if i == 1;
    plot(mean(topgrad_dat.GP60_L{1,i}', 2), 'Color', bluecol(6,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.GP60_R{1,i}', 2), 'Color', orgcol(6,:), 'LineWidth', 1.5)
    else
    plot(mean(topgrad_dat.GP60_L{1,i}', 2)+(i-1)*yoffset, 'Color', bluecol(6,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.GP60_R{1,i}', 2)+(i-1)*yoffset, 'Color', orgcol(6,:), 'LineWidth', 1.5)
    end
end

ylim([0 1.3*10^-11])

xline(101, '--k')
xticks([1:20:165])
xticklabels([])
xlim([41 161])

ylabel('Relative ERF amp');
set(gca, 'FontSize', txtsize)

ax = gca;
ax.XGrid = 'on';

plot([43 53 67 77 79 89 91 101], [0.9 0.9 0.61 0.61 0.38 0.38 0.1 0.1]*10^-11, '|k', 'MarkerSize', 22)

nexttile(3,1:2); hold on;

patch([112 112 132 132], [0 1.3*10^-11 1.3*10^-11 0], [0.2 0.2 0.2], 'FaceAlpha', 0.125, 'EdgeAlpha', 0)

plot(mean(topgrad_dat.PO60_R{1,5}', 2), 'Color', orgcol(6,:), 'LineWidth', 1.5)
plot(mean(topgrad_dat.PO60_L{1,5}', 2), 'Color', bluecol(6,:), 'LineWidth', 1.5)

%GPISI240
plot(mean(topgrad_dat.GP60_R{1,4}', 2), ':', 'Color', orgcol(6,:), 'LineWidth', 1.5)
plot(mean(topgrad_dat.GP60_L{1,4}', 2), ':', 'Color', bluecol(6,:), 'LineWidth', 1.5)

legend({'', 'Pulse only RIGHT', 'Pulse only LEFT', 'Gap+Pulse (isi 240ms) RIGHT', 'Gap+Pulse (isi 240ms) LEFT'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

ylim([0 1.3*10^-11])

xline(101, '--k')
xticks([1:20:165])
xticklabels([-500:100:320])
xlim([41 161])

ylabel('ERF amp (T/cm)');
xlabel('Time (ms relative pulse onset)')
set(gca, 'FontSize', txtsize)

ax = gca;
ax.XGrid = 'on';

%PO box L
nexttile(7); hold on;
boxchart(amp.PO60_L(:,:), 'BoxFaceColor', blue, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', blue, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'70', '75', '80', '85', '90', '95'})
ylabel('N1 ERF amp (T/cm)')
xlabel('Pulse level (dB)')
set(gca, 'FontSize', txtsize)

%PO box R
nexttile(10); hold on;
boxchart(amp.PO60_R(:,:), 'BoxFaceColor', orange, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', orange, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'70', '75', '80', '85', '90', '95'})
xlabel('Pulse level (dB)')
set(gca, 'FontSize', txtsize)

%GP Box L
nexttile(8); hold on;
boxchart(amp.GP60_L(:,:), 'BoxFaceColor', blue, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', blue, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'0', '60', '120', '240'})
ylabel('N1 ERF amp (T/cm)')
xlabel('ISI (ms)')
set(gca, 'FontSize', txtsize)

%GP box R
nexttile(11); hold on;
boxchart(amp.GP60_R(:,:), 'BoxFaceColor', orange, 'BoxFaceAlpha', 1, 'BoxLineColor', [1 1 1], 'MarkerColor', orange, 'MarkerSize', 3)
ylim([0 1.5*10^-11])
xticklabels({'0', '60', '120', '240'})
xlabel('ISI (ms)')
set(gca, 'FontSize', txtsize)

saveas(gcf, '../Analysis Output/Fig2_form9.svg')
close;
