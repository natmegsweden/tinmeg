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

[count60L, name60L] = groupcounts(L_topgrads.PO60chan);
[count60R, name60R] = groupcounts(R_topgrads.PO60chan);

[count70L, name70L] = groupcounts(L_topgrads.PO70chan);
[count70R, name70R] = groupcounts(R_topgrads.PO70chan);

%% Gather response from all conditions from SOIs identified for POXX_90

topgrad_dat = struct();

for ii = 1:numel(conditions); %For condition PO60 and GP60 (index [1,3] in var: 'conditions')

    for iii = 1:numel(cond.([conditions{ii} 'label']));
    
            %The condition and stimuli
            name = cond.([conditions{ii} 'label']){iii};
            disp(name); %helps to monitor what's going on
    
        for i = 1:numel(sub_date.ID)
            
            if ii == 1 || ii == 3 || ii == 5 && iii == 1 %if 60 dB carrier
                
                temp = tlk_sub_cmb.(conditions{ii}){i,iii};
                topgrad_dat.([conditions{ii} '_L']){1,iii}(i,:) = temp(L_topgrads.PO60chind{i},:); 

                temp = tlk_sub_cmb.(conditions{ii}){i,iii};
                topgrad_dat.([conditions{ii} '_R']){1,iii}(i,:) = temp(R_topgrads.PO60chind{i},:);

            elseif ii == 2 || ii == 4 || ii == 5 && iii == 2 %if 70 dB carrier
                
                temp = tlk_sub_cmb.(conditions{ii}){i,iii};
                topgrad_dat.([conditions{ii} '_L']){1,iii}(i,:) = temp(L_topgrads.PO70chind{i},:); 

                temp = tlk_sub_cmb.(conditions{ii}){i,iii};
                topgrad_dat.([conditions{ii} '_R']){1,iii}(i,:) = temp(R_topgrads.PO70chind{i},:);

            end
        %For subject
        end

    %For stims 
    end

%For conditions
end


%% Figures of top gradiometers

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
figure; hold on;
for i = 1:6 % Six stim (0-240ms ISI) in GP60 condition
    plot(mean(topgrad_dat.PO60_L{1,i}', 2), 'Color', bluecol(i,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.PO60_R{1,i}', 2), 'Color', orgcol(i,:), 'LineWidth', 1.5)
end

%legend({}, 'Location', 'northwest');

xlabel('Time (ms)')
ylabel('Top gradiometer ERF amplitude')

xline(100, '--k')

% xline(111, 'k')
% xline(131, 'k')

xticks([0:20:165])
xticklabels([-500:100:320])
xlim([40 160])

ylim([0 1.3*10^-11])

%saveas(gcf, '../Analysis Output/L_R_POtest.svg')


%GP60
yoffset = 0.25*10^-11;

figure; hold on;
for i = 1:4 % Six stim (70-95db Pulse) in PO60 condition
    if i == 1;
    plot(mean(topgrad_dat.GP60_L{1,i}', 2), 'Color', bluecol(6,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.GP60_R{1,i}', 2), 'Color', orgcol(6,:), 'LineWidth', 1.5)
    else
    plot(mean(topgrad_dat.GP60_L{1,i}', 2)+(i-1)*yoffset, 'Color', bluecol(6,:), 'LineWidth', 1.5)
    plot(mean(topgrad_dat.GP60_R{1,i}', 2)+(i-1)*yoffset, 'Color', orgcol(6,:), 'LineWidth', 1.5)
    end
end

ylim([0 1.3*10^-11])

xline(100, '--k')
xticks([0:20:165])
xticklabels([-500:100:320])
xlim([40 160])

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

%How many n of most common grad
mostgrad = max([count60L; count60R; count70L; count70R]);

clear bluecol orgcol
%Blue color gradient
bluecol(:,1) = linspace(185, 1, mostgrad) ./256; % R
bluecol(:,2) = linspace(202, 55, mostgrad) ./256; % G
bluecol(:,3) = linspace(254, 203, mostgrad) ./256; % B

%Orange color gradient
orgcol(:,1) = linspace(241, 252, mostgrad) ./256; % R
orgcol(:,2) = linspace(215, 96, mostgrad) ./256; % G
orgcol(:,3) = linspace(177, 10, mostgrad) ./256; % B

%Color left blue
for i = 1:numel(count60L)
    [M,I] = sort(count60L, 'descend')
    name60L{I(1)}
    count60L(I(1))
    
    %Find sensor label and put gradient according to n of sensor in colorspace
    sensind = find(ismember(sensors.label, name60L{I(i)}(1:7))); %Search for first gradiometer
    colors2(sensind,:) = bluecol(count60L(I(i)),:) %Color first gradiometer
    colors2(sensind-1,:) = bluecol(count60L(I(i)),:) %Color other magnetometer
    colors2(sensind+1,:) = bluecol(count60L(I(i)),:) %Color other gradiomter
end

%Color right orange
for i = 1:numel(count60R)
    [M,I] = sort(count60R, 'descend')
    name60R{I(1)}
    count60R(I(1))
    
    %Find sensor label and put gradient according to n of sensor in colorspace
    sensind = find(ismember(sensors.label, name60R{I(i)}(1:7))); %Search for first gradiometer
    colors2(sensind,:) = orgcol(count60R(I(i)),:) %Color first gradiometer
    colors2(sensind-1,:) = orgcol(count60R(I(i)),:) %Color other magnetometer
    colors2(sensind+1,:) = orgcol(count60R(I(i)),:) %Color other gradiomter
end

figure('Position', [400 200 1800 1000]); hold on
subplot(1,2,1)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.9);
%ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

subplot(1,2,2)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.9);
%ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([-100 25])

%Some suspect sensors
find(ismember(L_topgrads.PO60chan, 'MEG2042+2043'))
find(ismember(L_topgrads.PO60chan, 'MEG0122+0123'))
find(ismember(L_topgrads.PO60chan, 'MEG0342+0343'))

find(ismember(R_topgrads.PO60chan, 'MEG2032+2033'))
find(ismember(R_topgrads.PO60chan, 'MEG1212+1213'))

%inspect others
%topoplot
%PO vs GP
%quantify

%inspect sub 4
figure; hold on;
plot(tlk_sub_cmb.PO60{4,5}(8,:)) %favvisen (0242)
plot(tlk_sub_cmb.PO60{4,5}(12,:)) %irl
