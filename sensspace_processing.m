
% Run setup for common variables if not already loaded
if exist('sub_date', 'var') == 0; run Conditions_triggers.m; end;

%% Process timelockeds per condition, gather subject - and grand average

%Structure for subject timleockeds.avg
tlk_sub = struct();

%Structure for subject timleockeds.avg with combined planar gradiometers
tlk_sub_cmb = struct();

for ii = 1:length(conditions)

trig = cond.([conditions{ii} 'trig']);
nstim = numel(trig);

    for iii = 1:nstim

        for i = 1:numel(sub_date.ID)

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
        cfg.preproc.hpfilter = 'no';
        
        cfg.trials = tempdat.trialinfo == trig(iii);
        
        timelockeds = ft_timelockanalysis(cfg, tempdat);
        
        clear tempdat
        
        %save individual timelockeds
        save(['../mat_data/timelockeds/ID' sub_date.ID{i} '/' (cond.(([conditions{ii} 'label'])){iii}) '_tlks.mat'], 'timelockeds');
                
        cfg = [];
        timelockeds_cmb = ft_combineplanar(cfg, timelockeds);

        %save individual timelockeds with combined planar
        save(['../mat_data/timelockeds/ID' sub_date.ID{i} '/' (cond.(([conditions{ii} 'label'])){iii}) '_tlks_cmb.mat'], 'timelockeds_cmb');
        
        %Write averages to structs for easy access
        tlk_sub_cmb.(conditions{ii}){i, iii} = timelockeds_cmb.avg;
        tlk_sub.(conditions{ii}){i, iii} = timelockeds.avg;
        
        end
        
    end

end

clear i ii iii trig nstim

save('../mat_data/timelockeds/tlk_sub.mat', 'tlk_sub');
save('../mat_data/timelockeds/tlk_sub_cmb.mat', 'tlk_sub_cmb');

%% Find 8 combined gradiometers with biggest response in grand average of subjects in PO60/70 90-trials

%WIP WIP WIP

%How many combined grads of interest?
n_top_chan = 8;

%For 60dB Carrier
%Create empty cell array
max_grad = cell(1,n_top_chan);

load 

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



