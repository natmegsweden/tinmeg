%% Load conditions for subjects, average trials

%NB, fs_ds = 200 is hardcoded in for adapting TOI to gap in identifying amp and lat

all_cmb_avg = load('../mat_data/timelockeds/all_cmb_avg.mat');
all_cmb_avg = all_cmb_avg.all_cmb_avg;

gravg_cmb = load('../mat_data/timelockeds/grand_avg_cmb.mat');
gravg_cmb = gravg_cmb.gravg_cmb;

mean_sub = load(['../mat_data/timelockeds/mean_sub.mat']);
mean_sub = mean_sub.timelockeds_cmb;
mean_sub.time = round(mean_sub.time,3); %Round to avoid issues with floating point precision in plots

%Load sensor positions from one subjects data
sensors = load('../mat_data/preprocessing/ID0539/PO60_ds_clean.mat');
sensors = sensors.cleaned4mat.grad;

%find n = 8(?) top_chan in control condition (PO60/70 90)
%plot amplitude of topchans
%boxplot of top_chans
%helmet
%permutations test

%% Find 8 combined gradiometers with biggest response in grand average of subjects in PO60/70 90-trials

%How many combined grads of interest?
n_top_chan = 8;

%For 60dB Carrier
%Create empty cell array
max_grad = cell(1,n_top_chan);

%Overwrite mean_sub with grand average to maintain "labels"
mean_sub.avg = gravg_cmb.PO60{5}; %col 5 = PO6090

%find max of first 102 chan (grads) at mean of samples 116:131 (75-150ms)
[val, ind] = sort(mean(mean_sub.avg(1:102,116:131),2), 'descend');

%write top i cmb_grads to struct
for i = 1:n_top_chan
    max_grad{1,i} = mean_sub.label{ind(i)};
end

top_chan60 = unique(max_grad)';

%For 70dB Carrier
%Create empty cell array
max_grad = cell(1,n_top_chan);

%Overwrite mean_sub with grand average to maintain "labels"
mean_sub.avg = gravg_cmb.PO70{4}; %col 5 = PO7090

%find max of first 102 chan (grads) at mean of samples 116:131 (75-150ms)
[val, ind] = sort(mean(mean_sub.avg(1:102,116:131),2), 'descend');

%write top i cmb_grads to struct
for i = 1:n_top_chan
    max_grad{1,i} = mean_sub.label{ind(i)};
end

top_chan70 = unique(max_grad)';

%% Sensshape plot with highlighted chips

cfg = [];
cfg.output = 'scalp';

%Segment template MRI for head mesh
mri_segmented = ft_volumesegment(cfg, template_mri);

%Transform to neuromag coordsys - same as for sensors
mri_segmented = ft_convert_coordsys(mri_segmented, 'neuromag');

%Create mesh skull
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);

%Move scalp to not clip through sensors and be positioned nicely in helmet
mesh_scalp.pos(:,3) = mesh_scalp.pos(:,3) - 35; %Down
mesh_scalp.pos(:,2) = mesh_scalp.pos(:,2) - 15; %Back

%Plot higlighted sensors
%Write colour vector to row in mean_sub.label (+/-1 to row for all sensors on chip)
colors = ones(306, 3);

%Index based on the first gradiometer in top_chan
for i = 1:numel(top_chan60)
top_chan_temp{i} = top_chan60{i}(1:end-5);
end
top_chan_idx = ismember(sensors.label, top_chan_temp);
clear top_chan_temp

%match "colors" to index
for i = 1:numel(sensors.label)
    if top_chan_idx(i) == 1
    colors(i, :) = [0.85 0.325 0.098];
    colors(i+1, :) = [0.85 0.325 0.098];
    colors(i-1, :) = [0.85 0.325 0.098];
    end
end

[top_chan60 top_chan70]
warning('Note that sensors are the same for PO60 and PO70 - 90dB when n of sensors = 8');

figure('Position', [400 200 1800 1000]); hold on
subplot(1,2,1)
ft_plot_sens(sensors, 'facecolor', colors, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

subplot(1,2,2)
ft_plot_sens(sensors, 'facecolor', colors, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([-100 25])

%% Quantify SOI amplitude

%% Manuscript figures

