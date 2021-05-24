%% Read subjects and dates

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

disp(['Number of subjects in table is ' num2str(height(sub_date))])

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

%% Specify conditions and event triggers

run('Conditions_triggers.m');

%% Read or create log file of conditions

%create log for n of trials in raw data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_raw.csv', 'file');
    rawcondlog = readtable('../Analysis Output/n_cond_raw.csv', 'ReadVariableNames', false);
    rawcondlog = table2cell(rawcondlog);
else 
    %NB 22 columns hardcoded
    rawcondlogl = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

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