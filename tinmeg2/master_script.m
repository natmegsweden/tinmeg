
%% Specify conditions, event triggers and subject list
run('cond_trig_sub_tinmeg2.m');

meg_data_path = '/archive/20061_tinnitus/MEG/';

if ~exist('tlk_sub_cmb', 'var')
    tlk_sub_cmb = load('../mat_data/timelockeds/tinmeg2/tlk_sub_cmb.mat');
    tlk_sub_cmb = tlk_sub_cmb.tlk_sub_cmb;
end

%% Find tinmeg2-files for subject and create cell-array of file paths

% Create cell array for subjects filepaths
subpaths = cell(1);

for i = 1:1%height(sub_date);

% Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
fnames = find_files(subpath, {'tinmeg2', 'tsss'}, 'ds');

subpaths{i,1} = char(sub_date.ID{i}); %Include ID for tracking

    for fileindex = 1:length(fnames);
        subpaths{i,1+fileindex} = [subpath char(fnames(fileindex))]; % NB! n of files differ between rows, some subjects have empty columns
    end
clear fnames subpath i fileindex
end

%writetable(cell2table(subpaths), '../Analysis Output/included_filepaths_tinmeg2.csv') %Write log

%% create log for n of trials in raw data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_raw.csv', 'file');
    rawcondlog = readtable('../Analysis Output/n_cond_raw.csv', 'ReadVariableNames', false);
    rawcondlog = table2cell(rawcondlog);
else 
    %NB 22 columns hardcoded
    rawcondlog = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

%% Preprocess - loop for subject list in Full analysis
i = 2;
% specify i for now
fs_ds = 200;

%Check if ID is in TRIAL-log and determine row in log to write to
if find(strcmp(['ID' char(subpaths(i,1))], rawcondlog)) > 0;
    rawlogheight = find(strcmp(['ID' char(subpaths(i,1))], rawcondlog));
else
    %Find height of trial-log and +1 for new row
    rawlogheight = size(rawcondlog);
    rawlogheight = rawlogheight(1) + 1;

    %Write ID to new row, column 1
    rawcondlog{rawlogheight,1} = ['ID' char(subpaths(i,1))];
end

outdir = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/']

%Check if subject dir exist, create/define
if ~exist(outdir, 'file');
mkdir(outdir);
end

%For each condition
for ii = 1:length(conditions)

%Create filename and check if raw file for condition exist
fname = [conditions{ii} '_raw' '.mat'];

    if exist([outdir fname], 'file')
    warning(['Output ' fname ' exist for subject ' sub_date.ID{i}])
    continue
    end

% Define trials
cfg = [];

sub_data = subpaths(i,2:max(find(~cellfun(@isempty,subpaths(i,:)))));

% NB! cellfun for cfg.dataset defines 2:highest populated column
cfg.dataset             = sub_data;
cfg.trialdef.prestim    = 0.500;        % seconds before trigger
cfg.trialdef.poststim   = 0.500;        % seconds after trigger
cfg.trialdef.eventvalue = eval(['cond.' char(conditions(ii)) 'trig']); % :/
cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';

% On pre/post-stim, consider ISI240 as longest trial, i.e trial structure of
% baseline: 200 ms
% gap duration: 50 ms
% ISI: 240 ms, i.e, min of: 490 ms prestim

cfg = ft_definetrial(cfg);

%Remove trials from cfg.trl that have negative sample index for trial start
cfg.trl = cfg.trl(cfg.trl(:,1) >= 0,:);

%Remove trials from cfg.trl that have higher sample index than exist in file
cfg.trl = cfg.trl(cfg.trl(:,2) < max([cfg.event.sample]),:);

% preprocessing
cfg.demean     = 'no';
cfg.lpfilter   = 'no';
cfg.hpfilter   = 'no';
cfg.dftfilter  = 'no';
cfg.channel    = {'MEG', 'ECG', 'EOG'};

pre_tinmeg2 = ft_preprocessing(cfg);

save([outdir fname], 'pre_tinmeg2', '-v7.3')

%write n of trials log        
%How many stimtypes for cond and what trigger values
nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
trigs = eval(['cond.' char(conditions(ii)) 'trig']);

    %omg this gets ugly..
    %Write stim to rawcondlog cell-array
    for iii = 1:length(eval(['cond.' char(conditions(ii)) 'trig']))
    if strcmp(conditions(ii), 'PO60')
    rawcondlog(rawlogheight,iii+1) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'PO70')
    rawcondlog(rawlogheight,iii+7) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP60')
    rawcondlog(rawlogheight,iii+12) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP70')
    rawcondlog(rawlogheight,iii+16) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GO')
    rawcondlog(rawlogheight,iii+20) = num2cell(sum(res4mat.trialinfo(:) == trigs(iii)));
    end
    end

%Write log to csv
writetable(cell2table(rawcondlog), '../Analysis Output/n_cond_raw_tinmeg2.csv', 'WriteVariableNames', false) %Write log

%downsample and save
cfg = [];
cfg.resamplefs = fs_ds;

pre_tinmeg2_ds = ft_resampledata(cfg, pre_tinmeg2);

%NB new fname
fname = [char(conditions(ii)) '_ds' '.mat'];

save([outdir fname], 'pre_tinmeg2_ds')

%clear temp variables
clear pre_tinmeg2 pre_tinmeg2_ds

%end for conditions
end

%% Clean and log

%create log for n of trials in cleaned data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_clean.csv', 'file');
    cleancondlog = readtable('../Analysis Output/n_cond_clean.csv', 'ReadVariableNames', false);
    cleancondlog = table2cell(cleancondlog);
else 
    %NB 22 columns hardcoded
    cleancondlog = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

%% Loop over dowsampled data in ft_rejectvisual where no cleaned output

%Check if ID is in trial-log and determine row in log to write to
if find(strcmp(['ID' sub_date.ID{i}], cleancondlog)) > 0;
   logheight = find(strcmp(['ID' sub_date.ID{i}], cleancondlog));
else

%Find height of trial-log and +1 for new row
logheight = size(cleancondlog);
logheight = logheight(1) + 1;

%Write ID to new row, column 1
cleancondlog{logheight,1} = ['ID' sub_date.ID{i}];
end

%For each condition loaded
for ii = 1:length(conditions);
fname = [char(conditions(ii)) '_ds_clean' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

    %check if file exist
    if exist(fpath, 'file')
    warning([char(conditions(ii)) ' for subject: ID' sub_date.ID{i} ' exist'])
    continue
    end

%if not exist load downsampled (_ds) mat-file
fname = [char(conditions(ii)) '_ds' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

pre_tinmeg2_ds = load(fpath);
pre_tinmeg2_ds = pre_tinmeg2_ds.pre_tinmeg2_ds;

%ft_rejectvisual for condition
cfg = [];
cfg.method = 'summary';
cfg.keepchannel = 'yes';
cfg.channel = 'MEGMAG';
cfg.layout = 'neuromag306all.lay';

pre_tinmeg2_c = ft_rejectvisual(cfg,pre_tinmeg2_ds);

cfg.channel = 'MEGGRAD';

pre_tinmeg2_c = ft_rejectvisual(cfg, pre_tinmeg2_c);

%new filename and save
fname = [char(conditions(ii)) '_ds_clean' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

save(fpath, 'pre_tinmeg2_c');

%     %write n of trials log        
%     %How many stimtypes for cond and what trigger values
%     nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
%     trigs = eval(['cond.' char(conditions(ii)) 'trig']);
% 
%     %Write log for n of trials in condition after cleaning
%     for iii = 1:length(eval(['cond.' char(conditions(ii)) 'trig']))
%     if strcmp(conditions(ii), 'PO60')
%     cleancondlog(logheight,iii+1) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
%     elseif strcmp(conditions(ii), 'PO70')
%     cleancondlog(logheight,iii+7) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
%     elseif strcmp(conditions(ii), 'GP60')
%     cleancondlog(logheight,iii+12) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
%     elseif strcmp(conditions(ii), 'GP70')
%     cleancondlog(logheight,iii+16) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
%     elseif strcmp(conditions(ii), 'GO')
%     cleancondlog(logheight,iii+20) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
%     end
%     end
% 
% %Write log to csv
% writetable(cell2table(cleancondlog), '../Analysis Output/n_cond_clean.csv', 'WriteVariableNames', false) %Write log

%for conditions
end



%% ICA

for i = sub_date{i}
    run('tinmeg2'/ICA.m')
end

%% Create timelockeds

for i = sub_date{i}
    run('tinmeg2'/create_timelockeds.m')
end

%% Preliminary plot function

%sort highest response grads (1:102) at 75-150ms (116:131)
[val, ind] = sort(mean(tlk_sub_cmb.tin0_bkg0{1,3}(1:102,116:131),2), 'descend');

for j = [0,3,8] %tin

    for jjj = [0,3,8] %bkg
    
        name = ['tin' num2str(j) '_bkg' num2str(jjj)];
    
        figure('Position', [400 400 800 600]); tiledlayout(2,1, 'TileIndexing','columnmajor'); hold on;

        for jj = [3, 1]

            nexttile; hold on;
            
            %Eight top-grads
            for i = 1:8
                plot(tlk_sub_cmb.(name){1,jj}(ind(i), :), 'Color', [0 0 0 0.5]);
                temp(i,:) = tlk_sub_cmb.(name){1,jj}(ind(i), :);
            end
            plot(mean(temp), 'Color', [1 0 0], 'LineWidth', 2)
            clear temp
            
            if jj == 3;
                title([name '_PO'], 'Interpreter','none')
            elseif jj == 1;
                title([name '_GP'], 'Interpreter','none')
            end
            xlim([25 150])
            xline([101 101])
            xticks(0:25:201)
            xticklabels(-500:125:500)
        
        end
    
        saveas(gcf, ['../Analysis Output/' name '.svg'])
        close

    end

end

%% Mean plot only PO + GP

%WIP
%What subject from sub_date
sub = 1;

%sort highest response grads (1:102) at 75-150ms (116:131)
%[val, ind] = sort(mean(tlk_sub_cmb.tin0_bkg0{sub,3}(1:102,116:131),2), 'descend');

%TEMP for top_chan60
ind = find(ismember(timelockeds_cmb.label, top_chan60));

for j = [0,3,8] %tin

    for jjj = [0,3,8] %bkg
    
        name = ['tin' num2str(j) '_bkg' num2str(jjj)];
    
        figure('Position', [400 400 1100 400]); tiledlayout(2,1, 'TileIndexing','columnmajor'); hold on;

        for jj = [3] %Third cell in "tlk_sub_cmb" is PO

            nexttile; hold on;
            
            %Eight top-grads fo PO (i.e jj = 3)
            for i = 1:8
                %plot(tlk_sub_cmb.(name){1,jj}(ind(i), :), 'Color', [0 0 0 0.5]);
                temp(i,:) = tlk_sub_cmb.(name){sub,jj}(ind(i), :);
            end
            plot(mean(temp), 'LineWidth', 1.5)
            disp(['Max in PO avg ' name ': ' num2str(max(mean(temp)))])
            POmax=max(mean(temp(100:201)));
            clear temp

            %Eight top-grads for GPP (i.e jj = 1 or 3-1)
            for i = 1:8
                %plot(tlk_sub_cmb.(name){1,jj}(ind(i), :), 'Color', [0 0 0 0.5]);
                temp(i,:) = tlk_sub_cmb.(name){sub,jj-2}(ind(i), :);
            end
            plot(mean(temp), 'LineWidth', 1.5)
            disp(['Max in GP avg ' name ': ' num2str(max(mean(temp)))])
            GPmax=max(mean(temp(100:201)));
            clear temp
            
            title([name ' PO vs GP'], 'Interpreter','none')

            legend({"PO", "GP"}, 'Location', 'northwest', 'AutoUpdate', 'off');
            legend('boxoff');

            xlim([25 150])
            xline([100])
            xline([116 131])
            xticks(0:25:201)
            xticklabels(-500:125:500)

            disp(num2str(POmax-GPmax))
        
        end
        saveas(gcf, ['../Analysis Output/' name '_ica_avg_top60.svg'])
        close

    end

end

%% Sensshape plot with highlighted chips

%From TinMEG1
load('top_chan60.mat');

%load timelockeds
load('../mat_data/timelockeds/ID0905/PO_00_tlks.mat');
load('../mat_data/timelockeds/ID0905/PO_00_tlks_cmb.mat');

%Remove EOG/ECG channels to simplify plot sensors
cfg = [];
cfg.channel = 'MEG';
timelockeds = ft_selectdata(cfg,timelockeds);

%List sensors
sensors = timelockeds.grad;

[val, ind] = sort(mean(timelockeds_cmb.avg(1:102,116:131),2), 'descend');

topchan = timelockeds_cmb.label(ind(1:8));

%If not variable, load and segment mri and create skull mesh
if ~exist('mri_segmented', 'var')
    load standard_mri;
    template_mri = mri;
    
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

    clear mri
end

%Plot higlighted sensors
%Write colour vector to row in mean_sub.label (+/-1 to row for all sensors on chip)
colors2 = ones(306, 3);

%Index based on the first gradiometer in top_chan
%Sensshape needs idx for all 306 sensors which complicates this a bit.
for i = 1:numel(topchan)
top_chan_temp{i} = topchan{i}(1:end-5);
end
top_chan_idx = ismember(timelockeds.label, top_chan_temp);
clear top_chan_temp

%match "colors" to index
for i = 1:numel(timelockeds.label)
    if top_chan_idx(i) == 1
    colors2(i, :) = [0.85 0.325 0.098];
    colors2(i+1, :) = [0.85 0.325 0.098];
    colors2(i-1, :) = [0.85 0.325 0.098];
    end
end

figure('Position', [400 200 1800 1000], 'Renderer','painters'); hold on
subplot(1,2,1)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([100 25])

subplot(1,2,2)
ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 0.7);
ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [173/256 216/256 230/256]);
view([-100 25])

saveas(gcf, '../Analysis Output/pilot2_00.svg');
close;

%% EOG - create timelockeds EOG002 only

i = 1
destdirectory = ['../mat_data/timelockeds/' 'ID' sub_date.ID{i} '/EOG/'];

for ii = 1:length(conditions)

%NB reading in EOG data from file BEFORE ft_rejectvisual
fname_in = [conditions{ii} '_ds' '.mat'];
fpath_in = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname_in];

%loads as 'res4mat_ds'
res4mat_ds = load(fpath_in);
res4mat_ds = res4mat_ds.pre_tinmeg2_ds;

%Find triggers and stims in condition
nstim = length(eval(['cond.' conditions{ii} 'trig']));
trigs = eval(['cond.' conditions{ii} 'trig']);

triglabel = cond.([conditions{ii} 'label']);

%Extract EOG-timelockeds for all stim types
    for stim_index = 1:nstim

    trigger = trigs(stim_index);
    fname_out = [triglabel{stim_index} '_eog' '.mat'];

    cfg = [];

    cfg.channel = {'EOG002'};

    cfg.covariance = 'yes';
    cfg.covariancewindow = 'prestim';
    cfg.keeptrials = 'no'; %if yes, no avg in output variable "timelockeds"?
    cfg.preproc.demean = 'yes';

    if any(regexp(triglabel{stim_index}, 'PO_*')) || any(regexp(triglabel{stim_index}, 'GO_*'))
        %Baseline window: 150ms before pulse/gap onset in PO/GO trials
        cfg.preproc.baselinewindow = [-0.150 0];
    elseif any(regexp(triglabel{stim_index}, 'PPP_*')) || any(regexp(triglabel{stim_index}, 'GPP_*'))
        %Baseline window: 150ms before gap/pp onset in GP/PP-trials
        cfg.preproc.baselinewindow = [-0.440 -0.290];
    end

    cfg.preproc.lpfilter = 'yes';
    cfg.preproc.lpfreq = 70;
    cfg.trials = res4mat_ds.trialinfo == trigger;

    eog_timelockeds = ft_timelockanalysis(cfg, res4mat_ds);

    save([destdirectory fname_out], 'eog_timelockeds');

    clear eog_timelockeds

    %run again with cfg.keeptrials = yes
    cfg.keeptrials = 'yes';
    eog_timelockeds_all = ft_timelockanalysis(cfg, res4mat_ds);

    fname_out = [triglabel{stim_index} '_eog_all' '.mat'];

    save([destdirectory fname_out], 'eog_timelockeds_all');

    clear eog_timelockeds_all

    %end all stimtypes
    end

%end all conditions
end

%% load and plot EOG

i = 1;

tin = [0, 3, 8];
bkg = [0, 3, 8];

for ii = 1:numel(tin)

    for iii = 1:numel(bkg);
        
        name = ['POvsGP_EOG_tin' num2str(tin(ii)) 'bkg' num2str(bkg(iii))];

        eog_timelockeds = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/EOG/GPP_' num2str(tin(ii)) num2str(bkg(iii)) '_eog.mat']);
        GPP = eog_timelockeds.eog_timelockeds;

        eog_timelockeds = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/EOG/GPG_' num2str(tin(ii)) num2str(bkg(iii)) '_eog.mat']);
        GPG = eog_timelockeds.eog_timelockeds;
        
        figure('Position', [600 400 1100 400]); hold on
        plot(GPP.avg, 'LineWidth', 1.5)
        plot(GPG.avg, 'LineWidth', 1.5)

        title({'PO vs GP, EOG', ['tin' num2str(tin(ii)) 'bkg' num2str(bkg(iii))]})
        legend({"PO", "GP"}, 'Location', 'northwest', 'AutoUpdate', 'off');
        legend('boxoff');

        xlim([25 150])
        xline([100])
        xticks(0:25:201)
        xticklabels(-500:125:500)

        %saveas(gcf, ['../Analysis Output/tinmeg2_eog_test/' name '.svg']);
        close;

    end
end

%% TinMEG1 for subject

PO6090 = load('../mat_data/timelockeds/ID0539/PO60_90_tlks_cmb.mat');
PO6090 = PO6090.timelockeds_cmb;

GP240 = load('../mat_data/timelockeds/ID0539/GP60_i240_tlks_cmb.mat');
GP240 = GP240.timelockeds_cmb;

%sort highest response grads for PO - (1:102) at 75-150ms (116:131)
[val, ind] = sort(mean(PO6090.avg(1:102,116:131),2), 'descend');

figure; hold on;
%Eight top-grads fo PO (i.e jj = 3)
for i = 1:8
    %plot(tlk_sub_cmb.(name){1,jj}(ind(i), :), 'Color', [0 0 0 0.5]);
    temp(i,:) = PO6090.avg(ind(i), :);
end

plot(mean(temp), 'LineWidth', 1.5); clear temp

%Eight top-grads fo PO (i.e jj = 3)
for i = 1:8
    %plot(tlk_sub_cmb.(name){1,jj}(ind(i), :), 'Color', [0 0 0 0.5]);
    temp(i,:) = GP240.avg(ind(i), :);
end

plot(mean(temp), 'LineWidth', 1.5); clear temp




%% GO comparison




