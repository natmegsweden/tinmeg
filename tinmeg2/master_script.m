
%% Specify conditions, event triggers and subject list
run('cond_trig_sub_tinmeg2.m');

meg_data_path = ['/archive/20061_tinnitus/MEG/'];

% if ~exist('tlk_cmb_avg', 'var')
%     tlk_cmb_avg = load('../mat_data/timelockeds/tinmeg2/tlk_cmb_avg.mat');
%     tlk_cmb_avg = tlk_cmb_avg.tlk_cmb_avg;
% end

%% Find tinmeg2-files for subject and create cell-array of file paths

% Create cell array for subjects filepaths
subpaths = cell(1);

for i = 1:numel(sub_date.ID);

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
i = 3;
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
%writetable(cell2table(rawcondlog), '../Analysis Output/n_cond_raw_tinmeg2.csv', 'WriteVariableNames', false) %Write log

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



%% MR preparation, requires some manual checks and inputs

%  Require some manual input for fiducials and coordsys
for i = 1:numel(sub_date.ID);
    
    mri_path_in = ['../MRI/' 'NatMEG_' sub_date.MRI_ID{i}];
    mri_path_out = ['../mat_data/MRI_mat_tinmeg2/ID' sub_date.ID{i} '/'];
    
    %Check if realigned mri exist for subject
    if exist([mri_path_out 'mri_realigned.mat'], 'file')
        warning(['Output mri_realign already exist for subject: ' sub_date.ID{i}])
    continue
    
    elseif ~exist(mri_path_out, 'file')
    mkdir(mri_path_out);

    end
    run('tinmeg2/MR_prep.m');
    
end

clear('mri_coordsys', 'mri_realigned_1', 'mri_realigned_vol_2', 'mri_realigned_vol_3', 'sensshape', 'headshape', 'dcmfile', 'dicom_path', 'MEGfile', 'mri_path_in', 'mri_path_out', 'i');

%% MR preparation step 2

for i = 1:numel(sub_date.ID);
    
    %Note paths change from MR step 1 - in/out same here
    sub_mri_path = ['../mat_data/MRI_mat_tinmeg2/ID' sub_date.ID{i} '/'];
    fname = [sub_mri_path 'meg_headmodel.mat'];

    %Check if headmodel exist for subject
    if exist(fname, 'file')
    warning(['Output' fname ' exist for subject ' sub_date.ID{i}])
    continue
    end
    
    run('MR_prep2.m');
    
end

clear('template_grid', 'binary_brain', 'binary_scalp', 'binary_skull', 'headmodel_meg', 'mesh_brain', 'mri_segmented', 'mri_segmented_2', 'mri_realigned_vol_3', 'mri_resliced', 'cfg');


%% ICA

for i = 1:numel(sub_date.ID)
    if exist(['../mat_data/ICA/ID' sub_date.ID{i}], 'dir')
       warning(['ICA data exist for subject: ID' sub_date.ID{i}])
       continue
    else
        run('tinmeg2/ICA.m')
    end
end

%% Create timelockeds

for i = 1:numel(sub_date.ID)
    if exist(['../../mat_data/timelockeds/ID' sub_date.ID{i} '/'], 'dir')
       warning(['Timelockeds exist for subject: ID' sub_date.ID{i}])
       continue
    else
        run('tinmeg2/create_timelockeds.m')
    end
end

%% Frequency analysis

% Run ft_freqanalysis

%%


%% Two step grad selection: PO amp -> Inhib level

sub = 2;

%Tableau medium 10 palette
% palette = [173 139 201;
%     168 120 110;
%     114 158 206;
%     255 158 74;
%     237 102 93;
%     103 191 92;
%     237 151 202;
%     162 162 162;
%     205 204 93;
%     109 204 218]/256;

% palette = ones(10,3);
% palette(:,2) = linspace(0,0.8,10);
% palette(:,3) = linspace(0,0.8,10);

palette = viridis(8);

%Load subject PO struct if you need it
PO_00 = load(['../mat_data/timelockeds/ID' num2str(sub_date.ID{sub}) '/PO_00_tlks_cmb.mat']);
PO_00 = PO_00.timelockeds_cmb;

for j = 0 %[0,3,8] %tin

    for jjj = [0,3,8] %bkg
    
        name = ['tin' num2str(j) '_bkg' num2str(jjj)];
        figure('Position', [400 400 1100 800]); tiledlayout(2,1, 'TileIndexing','columnmajor'); hold on;

        %sort highest response grads (1:102) following pulse onset
        [val, ind] = sort(mean(tlk_cmb_avg.(name){sub,3}(1:102,101:end),2), 'descend');
        
        topchan = PO_00.label(ind(1:8));
        topchan

        for jj = [3] % 3 is PO in tlk_cmb_avg.xxx{subject, 3}

            nexttile; hold on;
            
            % Set new default color palette
            set(gca, 'ColorOrder', palette);

            %Eight top-grads amplitude plot
            plot(tlk_cmb_avg.(name){sub,jj}(ind(1:8), :)', 'LineWidth', 1.5);

            legend(PO_00.label(ind(1:8)), 'Box','off', 'AutoUpdate','off', 'Location','northwest', 'FontSize', 8);

            %mean amplitude
            %plot(mean(tlk_cmb_avg.(name){sub,jj}(ind(1:8), :), 1), 'LineWidth', 1.5);
            
            title({['Pulse Only']; [name ': n = 8, max mean amp (0-500ms) gradiometers']}, 'Interpreter','none')

%             if jj == 3;
%                 title([name '_PO'], 'Interpreter','none')
%             elseif jj == 1;
%                 title([name '_GP'], 'Interpreter','none')
%             end

            xlim([21 201])
            ylim([0*10^-11 1.5*10^-11])
            xline([101 101])
            xticks(1:25:201)
            xticklabels(-500:125:500)
            
            % Inhib plot - !! jj (i.e cond col) is static
            nexttile; hold on;
            
            set(gca, 'ColorOrder', palette);

            %inhib = tlk_cmb_avg.(name){sub,3}(ind(1:8), :) - tlk_cmb_avg.(name){sub,1}(ind(1:8), :); %PO-GP
            PO = tlk_cmb_avg.(name){sub,3}(ind(1:8), :);
            GP = tlk_cmb_avg.(name){sub,1}(ind(1:8), :);
            inhib =  ((PO-GP)./PO) .* 100; % 100.*((PO-GP)./PO)

            plot(inhib', 'LineWidth', 1.5)

            title({['Inhibition 100.*((PO-GP)./PO)']; [name ': same gradiometers']}, 'Interpreter','none')
            
            yline([0 0])
            xlim([21 201])
            %ylim([-6*10^-12 7.5*10^-12])
            ylim([0 100])
            xline([101 101])
            xticks(1:25:201)
            xticklabels(-500:125:500)
            xlabel('Time (ms)')
            ylabel('% Inhibition')

            clear inhib


        end
    
        saveas(gcf, ['../Analysis Output/' name '_' num2str(sub) '_ratio.svg'])
        close

        % senshape test (move up to previous block)

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
            cfg.numvertices = 500;
            
            mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);
            
            %Move scalp to not clip through sensors and be positioned nicely in helmet
            mesh_scalp.pos(:,3) = mesh_scalp.pos(:,3) - 35; %Down
            mesh_scalp.pos(:,2) = mesh_scalp.pos(:,2) - 15; %Back
        
            clear mri
        end
        
        %If not variable, load sensors and labels
        if ~exist('sensors', 'var')
            tlks = load('../mat_data/timelockeds/ID0916/PO_00_tlks.mat');
            tlks = tlks.timelockeds;
        
            sensors = tlks.grad;
            clear tlks;
        end
        
        %Plot higlighted sensors
        %Write colour vector to row in mean_sub.label (+/-1 to row for all sensors on chip)
        colors2 = ones(306, 3);
        
        %Index based on the first gradiometer in top_chan
        %Sensshape needs idx for all 306 sensors which complicates this a bit.
        for i = 1:numel(topchan)
        top_chan_temp{i} = topchan{i}(1:end-5);
        end
        top_chan_idx = ismember(sensors.label, top_chan_temp);
        %clear top_chan_temp
        
        %match "colors" to index
        palind = 1;
        for i = 1:8 % 8 highest grad amplitudes in condition from above
            %find first sensor from cmb among all sensors (now alphabetical)
            colrow = find(ismember(sensors.label, PO_00.label{ind(i)}(1:7)));
            colors2(colrow, :) = palette(palind,:);
            colors2(colrow+1, :) = palette(palind,:);
            colors2(colrow-1, :) = palette(palind,:);
            palind = palind+1; %advance to next color triplet in palette
        end
        
        figure('Position', [400 200 1800 1000], 'Renderer','painters'); hold on
        subplot(1,2,1)
        ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 1, 'label', 'off');
        ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [0.6758 0.8438 0.8984]);
        view([100 25])
        
        subplot(1,2,2)
        ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 1, 'label', 'off');
        ft_plot_mesh(ft_convert_units(mesh_scalp, 'cm'), 'edgecolor', [0.6758 0.8438 0.8984]);
        view([-100 25])
        
        %saveas(gcf, ['../Analysis Output/' name '_' num2str(sub) '_abs_head.svg'])
        close

    end

end


%% Plot sensors separately

%If not variable, load sensors and labels
if ~exist('sensors', 'var')
    tlks = load('../mat_data/timelockeds/ID0916/PO_00_tlks.mat');
    tlks = tlks.timelockeds;

    sensors = tlks.grad;
    clear tlks;
end

sub = 1;

%Load subject PO struct if you need it
PO_00 = load(['../mat_data/timelockeds/ID' num2str(sub_date.ID{sub}) '/PO_00_tlks_cmb.mat']);
PO_00 = PO_00.timelockeds_cmb;

for j = 0 %[0,3,8] %tin

    for jjj = [0,3,8] %bkg
    
        name = ['tin' num2str(j) '_bkg' num2str(jjj)];
        figure('Position', [316,39,1783,1182], 'Renderer','painters'); tiledlayout(8,2, 'TileSpacing','tight', 'Padding','tight'); hold on;

        %sort highest response grads (1:102) following pulse onset
        [val, ind] = sort(mean(tlk_cmb_avg.(name){sub,3}(1:102,101:end),2), 'descend');
        
        topchan = PO_00.label(ind(1:8));
        topchan

        for jj = [3] % 3 is PO in tlk_cmb_avg.xxx{subject, 3}, 1 is GPP
            
            %For sensor
            for k = 1:8
                nexttile; hold on;
    
                %Eight top-grads amplitude plot
                plot(tlk_cmb_avg.(name){sub,jj}(ind(k), :)', 'LineWidth', 1.5);
                plot(tlk_cmb_avg.(name){sub,jj-2}(ind(k), :)', 'LineWidth', 1.5);

                % PO - GP plot
                %plot(tlk_cmb_avg.(name){sub,jj}(ind(k), :) - tlk_cmb_avg.(name){sub,jj-2}(ind(k), :), 'LineWidth', 1.5);
    
                legend({'Pulse', 'Gap+Pulse'}, 'Box','off', 'AutoUpdate','off', 'Location','northwest', 'FontSize', 10);
    
                %mean amplitude
                %plot(mean(tlk_cmb_avg.(name){sub,jj}(ind(1:8), :), 1), 'LineWidth', 1.5);
                
                title({[num2str(k) ': ' topchan{k}]}, 'Interpreter','none')
    
                xlim([21 201])
                ylim([0*10^-11 1.5*10^-11])
                xline([101 101])
                xticks(1:25:201)
                xticklabels(-500:125:500)       

                %Plot higlighted sensors
                %Write colour vector to row in mean_sub.label (+/-1 to row for all sensors on chip)
                colors2 = ones(306, 3);
                
                %Index based on the first gradiometer in top_chan
                %Sensshape needs idx for all 306 sensors which complicates this a bit.
                for i = 1:numel(topchan)
                top_chan_temp{i} = topchan{i}(1:end-5);
                end
                top_chan_idx = ismember(sensors.label, top_chan_temp);
                %clear top_chan_temp
                
                %match "colors" to index
                %find first sensor from cmb among all sensors (now alphabetical)
                colrow = find(ismember(sensors.label, PO_00.label{ind(k)}(1:7)));
                colors2(colrow, :) = [1 0 0];
                colors2(colrow+1, :) = [1 0 0];
                colors2(colrow-1, :) = [1 0 0];
                
                %omg
                if sub == 2 && jjj == 0
                    sides = ['L', 'L', 'R', 'R', 'R', 'R', 'R', 'R'];
                elseif sub == 2 && jjj == 3
                    sides = ['R', 'L', 'L', 'R', 'R', 'L', 'R', 'L'];
                elseif sub == 2 && jjj == 8
                    sides = ['L', 'L', 'L', 'R', 'R', 'L', 'L', 'L'];
                elseif sub == 1 && jjj == 0
                    sides = ['R', 'L', 'R', 'L', 'L', 'R', 'R', 'R'];
                elseif sub == 1 && jjj == 3
                    sides = ['R', 'R', 'R', 'R', 'L', 'R', 'L', 'R'];
                elseif sub == 1 && jjj == 8
                    sides = ['R', 'R', 'L', 'L', 'R', 'R', 'R', 'R'];
                end

                side = sides(k);

                if side == 'R'
                    nexttile; hold on;
                    ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 1, 'label', 'off');
                    view([100 25])
                elseif side == 'L'
                    nexttile; hold on;
                    ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 1, 'label', 'off');
                    view([-100 25])
                end

            %For sensor
            end

        end
    
%     saveas(gcf, ['../Analysis Output/' name '_' num2str(sub) '_sep_sens.svg'])
%     close;
        
    %For bkg
    end

%For tin
end

%% Plot L/R sensor of interest separately and MANUAL patch TOI

%If not variable, load sensors and labels
if ~exist('sensors', 'var')
    tlks = load('../mat_data/timelockeds/ID0916/PO_00_tlks.mat');
    tlks = tlks.timelockeds;

    sensors = tlks.grad;
    clear tlks;
end

sub = 3;
warning(['SUBJECT ID: ' sub_date.ID{sub}])

%Load subject PO struct if you need it
PO_00 = load(['../mat_data/timelockeds/ID' num2str(sub_date.ID{sub}) '/PO_00_tlks_cmb.mat']);
PO_00 = PO_00.timelockeds_cmb;

%Define SOI, right first!
if sub_date.ID{sub} == '0916'
    soi = {'MEG1332+1333', 'MEG0242+0243'};
elseif sub_date.ID{sub} == '0905'
    soi = {'MEG1332+1333', 'MEG1612+1613'};
end

%This bit is to maintain order within find(ismember())
[tf,loc] = ismember(PO_00.label, soi);
[~,p] = sort(loc(tf));
soi_ind = find(tf);
soi_ind = soi_ind(p);
clear tf loc p

%tois = struct();

ylims = [0*10^-11 1.5*10^-11];

for j = 0 %[0,3,8] %tin

    for jjj = [0,3,8] %bkg
    
        name = ['tin' num2str(j) '_bkg' num2str(jjj)];
        figure('Position', [316,39,1783,1182], 'Renderer','painters'); tiledlayout(8,2, 'TileSpacing','tight', 'Padding','tight'); hold on;

        %sort highest response grads (1:102) following pulse onset
        ind = find(ismember(PO_00.label, soi));

        for jj = [3] % 3 is PO in tlk_cmb_avg.xxx{subject, 3}, 1 is GPP
            
            %For sensor
            for k = 1:2
                nexttile; hold on;
                
                %Specify signal vectors for PO and GPP
                povec = tlk_cmb_avg.(name){sub,jj}(soi_ind(k), :);
                gpvec = tlk_cmb_avg.(name){sub,jj-2}(soi_ind(k), :);
                
                %specify tois manually
                if jjj == 0 && ismember(sub_date.ID(sub),'0916') && ismember(soi(k),'MEG1332+1333');
                    patch1on = 111;   patch1off = 127;
                    patch2on = 127;   patch2off = 140;
                elseif jjj == 0 && ismember(sub_date.ID(sub),'0916') && ismember(soi(k),'MEG0242+0243');
                    patch1on = 115;   patch1off = 128;
                    patch2on = 128;   patch2off = 146;
                elseif jjj == 3 && ismember(sub_date.ID(sub),'0916') && ismember(soi(k),'MEG1332+1333');
                    patch1on = 113;   patch1off = 127;
                    patch2on = 127;   patch2off = 137;
                elseif jjj == 3 && ismember(sub_date.ID(sub),'0916') && ismember(soi(k),'MEG0242+0243');
                    patch1on = 127;   patch1off = 141;
                    patch2on = 116;   patch2off = 127;
                elseif jjj == 8 && ismember(sub_date.ID(sub),'0916') && ismember(soi(k),'MEG1332+1333');
                    patch1on = 113;   patch1off = 127;
                    patch2on = 127;   patch2off = 136;
                elseif jjj == 8 && ismember(sub_date.ID(sub),'0916') && ismember(soi(k),'MEG0242+0243');
                    patch1on = 126;   patch1off = 144;
                    patch2on = 116;   patch2off = 126;
                elseif jjj == 0 && ismember(sub_date.ID(sub),'0905') && ismember(soi(k),'MEG1332+1333'); %%
                    patch1on = 114;   patch1off = 131;
                    patch2on = 131;   patch2off = 142;
                elseif jjj == 0 && ismember(sub_date.ID(sub),'0905') && ismember(soi(k),'MEG1612+1613');
                    patch1on = 115;   patch1off = 130;
                    patch2on = 130;   patch2off = 146;
                elseif jjj == 3 && ismember(sub_date.ID(sub),'0905') && ismember(soi(k),'MEG1332+1333');
                    patch1on = 114;   patch1off = 129;
                    patch2on = 129;   patch2off = 143;
                elseif jjj == 3 && ismember(sub_date.ID(sub),'0905') && ismember(soi(k),'MEG1612+1613');
                    patch1on = 127;   patch1off = 146;
                    patch2on = 116;   patch2off = 127;
                elseif jjj == 8 && ismember(sub_date.ID(sub),'0905') && ismember(soi(k),'MEG1332+1333');
                    patch1on = 116;   patch1off = 128;
                    patch2on = 128;   patch2off = 145;
                elseif jjj == 8 && ismember(sub_date.ID(sub),'0905') && ismember(soi(k),'MEG1612+1613');
                    patch1on = 127;   patch1off = 149;
                    patch2on = 117;   patch2off = 127;
                end
                
%                 %save right sensor tois to struct
%                 if ismember(soi(k),'MEG1332+1333')
%                     tois.([name '_' num2str(sub)]){1,1} = patch1on;
%                     tois.([name '_' num2str(sub)]){1,2} = patch1off;
%                     tois.([name '_' num2str(sub)]){2,1} = patch2on;
%                     tois.([name '_' num2str(sub)]){2,2} = patch2off;
%                 end
%     
%                 patch([patch1on patch1on patch1off patch1off], [ylims(1) ylims(2) ylims(2) ylims(1)], [252, 186, 18]/256, 'FaceAlpha', 0.7, 'EdgeAlpha', 0.3)
%                 patch([patch2on patch2on patch2off patch2off], [ylims(1) ylims(2) ylims(2) ylims(1)], [252, 186, 18]/256, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.3)
%                 %patch([patch3on patch3on patch3off patch3off], [ylims(1) ylims(2) ylims(2) ylims(1)], [253, 220, 34]/256, 'FaceAlpha', 0.75, 'EdgeAlpha', 0)
                
                %Plot responses
                plot(povec', 'Color', [35, 106, 185]/256, 'LineWidth', 1.5);
                plot(gpvec', 'Color', [214, 26, 70]/256, 'LineWidth', 1.5);

                legend({'TOI 1', 'TOI 2', 'Pulse', 'Gap+Pulse'}, 'Box','off', 'AutoUpdate','off', 'Location','northwest', 'FontSize', 10);

                %mean amplitude
                %plot(mean(tlk_cmb_avg.(name){sub,jj}(ind(1:8), :), 1), 'LineWidth', 1.5);
                
                title({[num2str(k) ': ' soi{k}]}, 'Interpreter','none');

                xlim([21 201])
                ylim([ylims(1) ylims(2)])
                xline([101 101])
                xticks(1:25:201)
                xticklabels(-500:125:500)

                %Plot higlighted sensors
                %Write colour vector to row in mean_sub.label (+/-1 to row for all sensors on chip)
                colors2 = ones(306, 3);
                
                %Index based on the first gradiometer in top_chan
                %Sensshape needs idx for all 306 sensors which complicates this a bit.
                for i = 1:numel(soi)
                top_chan_temp{i} = soi{i}(1:end-5);
                end
                top_chan_idx = ismember(sensors.label, top_chan_temp);
                %clear top_chan_temp
                
                %match "colors" to index
                %find first sensor from cmb among all sensors (now alphabetical)
                colrow = find(ismember(sensors.label, PO_00.label{soi_ind(k)}(1:7)));
                colors2(colrow, :) = [1 0 0];
                colors2(colrow+1, :) = [1 0 0];
                colors2(colrow-1, :) = [1 0 0];
                
                sides = ['L', 'R'];

                side = sides(k);

                if side == 'L'
                    nexttile; hold on;
                    ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 1, 'label', 'off');
                    view([100 25])
                elseif side == 'R'
                    nexttile; hold on;
                    ft_plot_sens(sensors, 'facecolor', colors2, 'facealpha', 1, 'label', 'off');
                    view([-100 25])
                end

            %For sensor
            end

        end
    
     saveas(gcf, ['../Analysis Output/' name '_' num2str(sub) '_sep_soi.svg'])
     close;
        
    %For bkg
    end

%For tin
end

%% Gather 1L & 1R sensor per subject for PO and GPP

L_R_soi = struct();

%For subject
for i = 1:numel(sub_date.ID);

    %Define SOI, right first!
    if sub_date.ID{i} == '0916'
        soi = {'MEG1332+1333', 'MEG0242+0243'};
    elseif sub_date.ID{i} == '0905'
        soi = {'MEG1332+1333', 'MEG1612+1613'};
    end

    %For bkg
    for ii = [0,3,8];
        
        %tin
        for iii = [0,3,8];

            fnamein = ['tin' num2str(iii) '_bkg' num2str(ii)];

            PO = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/PO_' num2str(iii) num2str(ii) '_tlks_cmb.mat']);
            PO = PO.timelockeds_cmb;

            GP = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/GPP_' num2str(iii) num2str(ii) '_tlks_cmb.mat']);
            GP = GP.timelockeds_cmb;
    
            %Keep only channels of interest
            cfg = [];
            cfg.channel = soi;
            PO = ft_selectdata(cfg, PO);
            GP = ft_selectdata(cfg, GP);

            L_R_soi.(['R_PO_' num2str(iii) num2str(ii)])(i,:) = PO.avg(find(ismember(PO.label, soi{1})),:);
            L_R_soi.(['L_PO_' num2str(iii) num2str(ii)])(i,:) = PO.avg(find(ismember(PO.label, soi{2})),:);
            L_R_soi.(['R_GP_' num2str(iii) num2str(ii)])(i,:) = GP.avg(find(ismember(GP.label, soi{1})),:);
            L_R_soi.(['L_GP_' num2str(iii) num2str(ii)])(i,:) = GP.avg(find(ismember(GP.label, soi{2})),:);

        %for tin
        end

    %for bkg    
    end

%for subject    
end

save('../mat_data/timelockeds/tinmeg2/L_R_soi.mat', 'L_R_soi')

%% Plot 1L & 1R sensor per subject for PO and GPP

if ~exist('L_R_soi', 'var')
    load('../mat_data/timelockeds/tinmeg2/L_R_soi.mat');
end

hemis = ['L', 'R'];

%For subject
for j = 1:numel(sub_date.ID)
    
    for i = 1:numel(hemis)
    
        figure('Position', [494,6,650,855]); tiledlayout(3,1,"TileSpacing","tight"); hold on;

%         %For tin
%         for ii = [0,3,8];
    
            %For bkg
            for iii = [0,3,8];
                
                %3x3
%                 PO = L_R_soi.([hemis(i) '_PO_' num2str(ii) num2str(iii)])(j,:);
%                 GP = L_R_soi.([hemis(i) '_GP_' num2str(ii) num2str(iii)])(j,:);
%     
%                 nexttile; hold on;
%                 plot(PO');
%                 plot(GP');
% 
%                 title([sub_date.ID{j} ' ' hemis(i) ' tin' num2str(ii) ' bkg' num2str(iii)]);

                PO_tin0 = L_R_soi.([hemis(i) '_PO_' num2str(0) num2str(iii)])(j,:);
                GP_tin0 = L_R_soi.([hemis(i) '_GP_' num2str(0) num2str(iii)])(j,:);
                tin0_in = PO_tin0-GP_tin0;

                PO_tin3 = L_R_soi.([hemis(i) '_PO_' num2str(3) num2str(iii)])(j,:);
                GP_tin3 = L_R_soi.([hemis(i) '_GP_' num2str(3) num2str(iii)])(j,:);
                tin3_in = PO_tin3-GP_tin3;

                PO_tin8 = L_R_soi.([hemis(i) '_PO_' num2str(8) num2str(iii)])(j,:);
                GP_tin8 = L_R_soi.([hemis(i) '_GP_' num2str(8) num2str(iii)])(j,:);
                tin8_in = PO_tin8-GP_tin8;

                nexttile; hold on;
                plot(tin0_in')
                plot(tin3_in')
                plot(tin8_in')

                title([sub_date.ID{j} ' ' hemis(i) ' bkg' num2str(iii)]);
                legend({'Tin 0', 'Tin 3', 'Tin 8'}, 'AutoUpdate', 'off', 'Location','southwest')

                xlim([21 201])
                ylim([-5*10^-12 7*10^-12])
                xline([101 101])
                yline([0 0])
                xticks(1:25:201)
                xticklabels(-500:125:500) 

            %For bkg
            end
    
%         %For tin
%         end

    saveas(gcf, ['../Analysis Output/' sub_date.ID{j} '_' hemis(i) '.svg']);
    close;

    %For hemis
    end

%For subject
end

%% Plot SOI for GO

%If not variable, load sensors and labels
if ~exist('sensors', 'var')
    tlks = load('../mat_data/timelockeds/ID0916/PO_00_tlks.mat');
    tlks = tlks.timelockeds;

    sensors = tlks.grad;
    clear tlks;
end

%for sub i !!
sub = 2;

i = sub;

%Load subject PO struct if you need it
% PO_00 = load(['../mat_data/timelockeds/ID' num2str(sub_date.ID{sub}) '/PO_00_tlks_cmb.mat']);
% PO_00 = PO_00.timelockeds_cmb;

%Define SOI, right first!
if sub_date.ID{i} == '0916'
    soi = {'MEG1332+1333', 'MEG0242+0243'};
elseif sub_date.ID{i} == '0905'
    soi = {'MEG1332+1333', 'MEG1612+1613'}; 
end

%This bit is to maintain order within find(ismember())
[tf,loc] = ismember(PO_00.label, soi);
[~,p] = sort(loc(tf));
soi_ind = find(tf);
soi_ind = soi_ind(p);
clear tf loc p

for j = 0 %[0,3,8] %tin

    figure('Position', [316,39,1000,1182], 'Renderer','painters'); tiledlayout(3,1, 'TileSpacing','tight', 'Padding','tight'); hold on;

    for jjj = [0,3,8] %bkg
    
        name = ['tin' num2str(j) '_bkg' num2str(jjj)];

        nexttile; hold on;

        for jj = [4] % 4 is condition GO in tlk_cmb_avg.xxx{subject, condition}
            
            %For sensor
            for k = soi_ind;

                %GO from soi_grad
                plot(tlk_cmb_avg.(name){sub,jj}(k, :)', 'LineWidth', 1.5);

                % PO - GP plot
                %plot(tlk_cmb_avg.(name){sub,jj}(ind(k), :) - tlk_cmb_avg.(name){sub,jj-2}(ind(k), :), 'LineWidth', 1.5);
    
                legend({[soi{1} ' (R)'], [soi{2} ' (L)']}, 'Box','off', 'AutoUpdate','off', 'Location','northwest', 'FontSize', 10);
    
                %mean amplitude
                %plot(mean(tlk_cmb_avg.(name){sub,jj}(ind(1:8), :), 1), 'LineWidth', 1.5);
                
                title({['SOI: GO, bgk: ', num2str(jjj)]}, 'Interpreter','none')
    
                xlim([21 201])
                ylim([0*10^-11 1.5*10^-11])
                xline([101 101])
                xticks(1:25:201)
                xticklabels(-500:125:500)       
                
                %match "colors" to index

            %For sensor
            end

        end
        
    %For bkg
    end

%For tin
end


    saveas(gcf, ['../Analysis Output/' name '_' num2str(sub) '_soi_GO.svg'])
    close;


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




