%% Read subjects and dates

%To do: 
%implement "run for who" to re-run single subject/testing purposes - stop
%if invalid ID

% Run for who?
runwho = input('Who are we analysing today? specify all (default) or specific ID \n', 's');
if isempty(runwho)
    runwho = 'all'
end

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%MRI data path
mri_data_path = '../MRI/';

%Sourcemodel template
load('/../../fieldtrip-20210311/template/sourcemodel/standard_sourcemodel3d6mm');
template_grid = sourcemodel;
clear sourcemodel

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

%% Read or create log file of conditions and artifacts

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

clear('i', 'ii', 'iii', 'logheight', 'fs_ds');

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

clear('i', 'ii', 'iii', 'logheight');

%% Timelockedanalysis



%% MR step 1
%  Require some manual input for fiducials and coordsys

for i = 1:length(sub_date.ID);
    
    sub_mri_path = ['../MRI/' 'NatMEG_' char(sub_date{i,1})];
    fname = ['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_realigned.mat'];

    %Check if realigned mri exist for subject
    if exist(fname, 'file')
    warning(['Output' fname ' exist for subject ' char(sub_date{i,1})])
    continue
    end
    
    run('D_MR_prep.m');
    
end

clear('');

%% MR Step 2
%  Re-slice and segment, time-conusming (5 min per subject) not requiring manual inputs

for i = 1:length(sub_date.ID);
    
    sub_mri_path = ['../MRI/' 'NatMEG_' char(sub_date{i,1})];
    fname = ['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_MEG_headmodel.mat'];

    %Check if headmodel exist for subject
    if exist(fname, 'file')
    warning(['Output' fname ' exist for subject ' char(sub_date{i,1})])
    continue
    end
    
    run('D_MR_prep2.m');
    
end

%% MR Step 3
%  Loop through plots for ALL[!] subjects to check for errors

for i = 1:length(sub_date.ID);
    
    load(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_headshape']);
    load(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_sensshape']);
    load(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_MEG_headmodel']);
    
    %Final plot - aligned MEG   
    fig = figure;
    hold on;
    ft_plot_sens(sensshape)
    ft_plot_headshape(headshape)
    ft_plot_headmodel(headmodel_meg)
    ft_plot_axes([], 'unit', 'cm');
    
    %Pause loop until figure is closed
    uiwait(fig);
    
end

%% Source reconstruction Beamformer


