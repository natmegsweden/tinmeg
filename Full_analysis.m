%% Read subjects and dates

%To do: 
%implement "run for who" to re-run single subject/testing purposes - stop
%if invalid ID

% Run for who?
% runwho = input('Who are we analysing today? specify all (default) or specific ID \n', 's');
% if isempty(runwho)
%     runwho = 'all'
% end

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%MRI data path
mri_data_path = '../MRI/';

%Sourcemodel template
load('/../../fieldtrip-20210311/template/sourcemodel/standard_sourcemodel3d6mm');
template_grid = sourcemodel;
template_grid = ft_convert_units(template_grid, 'mm');
clear sourcemodel;

%MRI template (https://www.fieldtriptoolbox.org/template/headmodel/)
load standard_mri;
template_mri = mri;
clear mri;

%Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

disp(['Number of subjects in table is ' num2str(height(sub_date))])

%Specify conditions and event triggers

run('Conditions_triggers.m');

%% Find relevant files for subject and create cell-array of file paths

% Create cell array for subjects filepaths
subpaths = cell(1);

for i = 1:height(sub_date);

% Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
fnames = find_files(subpath, {'tinmeg1', 'tsss'}, 'ds');

subpaths{i,1} = char(sub_date.ID{i}); %Include ID for tracking

    for fileindex = 1:length(fnames);
        subpaths{i,1+fileindex} = [subpath char(fnames(fileindex))]; % NB! n of files differ between rows, some subjects have empty columns
    end
clear fnames subpath i fileindex
end

writetable(cell2table(subpaths), '../Analysis Output/included_filepaths.csv') %Write log

%% Read or create log file of conditions and trials for preprocessing

%create log for n of trials in raw data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_raw.csv', 'file');
    rawcondlog = readtable('../Analysis Output/n_cond_raw.csv', 'ReadVariableNames', false);
    rawcondlog = table2cell(rawcondlog);
else 
    %NB 22 columns hardcoded
    rawcondlog = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

%create log for n of artifacts rejected
%first row are labels, first column IDs
% if exist('../Analysis Output/n_cond_artifact.csv', 'file');
%     artcondlog = readtable('../Analysis Output/n_cond_artifact.csv', 'ReadVariableNames', false);
%     artcondlog = table2cell(artcondlog);
% else 
%     %NB columns hardcoded
%     artcondlog = ['ID', strcat(conditions, '_jmp'), strcat(conditions, '_mus'), strcat(conditions, '_eog')];
% end

%%  Loop over A_Preprocess.m for subjects without output files
%   Also downsamples and saves with '_ds' in filename

%specify downsampled freq(Hz)
fs_ds = 200;

for i = 1:length(sub_date.ID);
    
    run('A_Preprocess.m');
    
end

clear('i', 'ii', 'iii', 'logheight', 'fs_ds', 'trigs', 'rawcondlog', 'fname', 'nstim', 'outdir', 'rawlogheight');

%% Loop over B_Clean.m for subjects and condititions with no cleaned output

%create log for n of trials in cleaned data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_clean.csv', 'file');
    cleancondlog = readtable('../Analysis Output/n_cond_clean.csv', 'ReadVariableNames', false);
    cleancondlog = table2cell(cleancondlog);
else 
    %NB 22 columns hardcoded
    cleancondlog = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

for i = 1:length(sub_date.ID);
    
    run('B_Clean.m');
    
end

clear('i', 'ii', 'iii', 'logheight', 'cleancondlog', 'cleaned4mat', 'fname', 'fpath', 'nstim', 'res4mat_ds', 'trigs');

%% Process timelockeds of EOG002

for i = 1:length(sub_date.ID)
    
    destdirectory = ['../mat_data/timelockeds/' 'ID' sub_date.ID{i} '/EOG/'];

    %Check if output folder exist and skip subject or create folder
    if exist(destdirectory, 'file');
       warning(['EOG output folder already exist for subject ' sub_date.ID{i} ' - skipping..']);
        continue
    elseif ~exist(destdirectory, 'file');
        mkdir(destdirectory);
    end
    
    run('C_process_EOG.m');
    
end

%% ICA

%Warning: Modified ft_artefact_ecg to assume 'y' for inspection 2021-09-17
% #205, #291, #306

for i = 1:length(sub_date.ID);
    
    outdir = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    
    if exist(outdir, 'file');
    warning(['Output folder exist for subject ' sub_date.ID{i} ' - skipping..']);
    continue
    end
    
    run('C2_ICA.m');
    
end


%% Timelockedanalysis




%% MR step 1
%  Require some manual input for fiducials and coordsys

for i = 1:4%length(sub_date.ID);
    
    mri_path_in = ['../MRI/' 'NatMEG_' sub_date.ID{i}];
    mri_path_out = ['../mat_data/MRI_mat/ID' sub_date.ID{i} '/'];
    
    %Check if realigned mri exist for subject
    if exist([mri_path_out 'mri_realigned.mat'], 'file')
    warning(['Output mri_realignet exist for subject: ' char(sub_date{i,1})])
    continue
    
    elseif ~exist(mri_path_out, 'file')
    mkdir(mri_path_out);

    end
    
    run('D_MR_prep.m');
    
end

clear('mri_coordsys', 'mri_realigned_1', 'mri_realigned_vol_2', 'mri_realigned_vol_3', 'submri', 'sensshape', 'headshape', 'dcmfile', 'dicom_path', 'MEGfile', 'mri_path_in', 'mri_path_out', 'i');

%% MR Step 2
%  Re-slice and segment, time-conusming (5 min per subject) not requiring manual inputs

for i = 1:4%length(sub_date.ID);
    
    %Note paths change from MR step 1 - in/out same here
    sub_mri_path = ['../mat_data/MRI_mat/ID' sub_date.ID{i} '/'];
    
    fname = [sub_mri_path 'meg_headmodel.mat'];

    %Check if headmodel exist for subject
    if exist(fname, 'file')
    warning(['Output' fname ' exist for subject ' sub_date.ID{i}])
    continue
    end
    
    run('D_MR_prep2.m');
    
end

clear('binary_brain', 'binary_scalp', 'binary_skull', 'headmodel_meg', 'mesh_brain', 'mri_segmented', 'mri_segmented_2', 'mri_realigned_vol_3', 'mri_resliced', 'cfg');

%% MR Step 3
%  Loop through plots for ALL[!] subjects to check for error

%Subject sourcemodel (basedonmni) and headmodel
for i = 1:4%length(sub_date.ID);
    meg_inpath = ['../mat_data/ICA/' 'ID' sub_date.ID{i} '/'];
    mri_inpath = ['../mat_data/MRI_mat/ID' sub_date.ID{i} '/'];
    outdir = ['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/'];
    
    %headmodel_meg
    headmodel_meg = load([mri_inpath 'meg_headmodel.mat']); 
    headmodel_meg = headmodel_meg.headmodel_meg;
    
    %Load subject sourcemodel in template grid format (based on MNI)
    subject_grid = load([mri_inpath 'subject_grid.mat']);
    subject_grid = subject_grid.subject_grid;
    
    
    %Final plot - aligned MEG   
    fig = figure('Position', [400 300 1800 800], 'Name', ['SUBJECT: ' sub_date.ID{i}], 'NumberTitle', 'off');
    hold on;
    
    subplot(1,2,1)
    ft_plot_headmodel(headmodel_meg, 'edgecolor', 'none', 'facealpha', 0.4,'facecolor', 'b');
    ft_plot_mesh(subject_grid.pos(subject_grid.inside,:));
    x = gca;
    x.CameraPosition = [-1100 -1300 400];
    
    subplot(1,2,2)
    ft_plot_headmodel(headmodel_meg, 'edgecolor', 'none', 'facealpha', 0.4,'facecolor', 'b');
    ft_plot_mesh(subject_grid.pos(subject_grid.inside,:));
    x = gca;
    x.CameraPosition = [1100 1300 200];
    
    %Pause loop until figure is closed
    uiwait(fig);
    
end

%Headshape and headmodel in sensors
for i = 1:4%length(sub_date.ID);
    
    sub_mri_path = ['../mat_data/MRI_mat/ID' sub_date.ID{i} '/'];
    
    %ft_read_headshape(MEGfile);
    load([sub_mri_path 'headshape']);
    
    %ft_read_sens(MEGfile);
    load([sub_mri_path 'sensshape']);
    
    %ft_prepare_headmodel(cfg, mesh_brain);
    load([sub_mri_path 'meg_headmodel']);
    
    %Final plot - aligned MEG
    fig = figure('Position', [400 300 1800 800], 'Name', ['SUBJECT: ' sub_date.ID{i}], 'NumberTitle', 'off'); %Postion: [Left Bottom Width Height]
    hold on;
    
    subplot(1,2,1)
    ft_plot_sens(sensshape)    
    ft_plot_headshape(headshape)
    ft_plot_headmodel(headmodel_meg)
    ft_plot_axes([], 'unit', 'cm');
    x = gca;
    x.CameraPosition = [240 140 140];
    
    subplot(1,2,2)
    ft_plot_sens(sensshape)
    ft_plot_headshape(headshape)
    ft_plot_headmodel(headmodel_meg)
    ft_plot_axes([], 'unit', 'cm');    
    
    x = gca;
    x.CameraPosition = [-240 -140 70];
    
    %Pause loop until figure is closed
    uiwait(fig);
    
end

%% Source reconstruction Beamformer


