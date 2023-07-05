%% Read paths, variables and templates

%Add functions to path
addpath('../functions/')

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%MRI data path
mri_data_path = '/archive/20061_tinnitus/MRI/';

%Sourcemodel template
load('../fieldtrip-20220711/template/sourcemodel/standard_sourcemodel3d6mm.mat')
template_grid = sourcemodel;
template_grid = ft_convert_units(template_grid, 'mm');
clear sourcemodel;

%MRI template (https://www.fieldtriptoolbox.org/template/headmodel/)
load standard_mri;
template_mri = mri;
clear mri;

%Read subject list
sub_date = readtable('../input_vars/sub_date.txt', 'Format', '%s%s%s%s%s');
disp(['Number of subjects in table is ' num2str(height(sub_date))])

%Read MEG triggers and conditions
run('../input_vars/triggers.m')

%Downsample rate
fs_ds = 200;

%% Find raw data files to include
run('find_fif.m')

%% Preprocess data to epochs

%Load subpaths if not a variable
if exist('subpaths','var') == 0
subpaths = readtable('../analysis_output/logs/included_filepaths.csv');
else

    %Run preprocessing for each subject
    for i = 1:numel(sub_date.ID)
    
        if sub_date.Exp{i} == 'tinmeg1'
           temp_cond = cond.tinmeg1.all;
           temp_stim = cond.tinmeg1;
        elseif sub_date.Exp{i} == 'tinmeg2'
           temp_cond = cond.tinmeg2.all;
           temp_stim = cond.tinmeg2;
        elseif sub_date.Exp{i} == 'tinmeg3'
           temp_cond = cond.tinmeg3.all;
           temp_stim = cond.tinmeg3;
        end

        run('preprocess.m')
    
    end

end

%% Loop over downsampled data and remove high variance trials/channels manually

%Run preprocessing for each subject
for i = 1:numel(sub_date.ID)

    if sub_date.Exp{i} == 'tinmeg1'
       temp_cond = cond.tinmeg1.all;
       temp_stim = cond.tinmeg1;
    elseif sub_date.Exp{i} == 'tinmeg2'
       temp_cond = cond.tinmeg2.all;
       temp_stim = cond.tinmeg2;
    elseif sub_date.Exp{i} == 'tinmeg3'
       temp_cond = cond.tinmeg3.all;
       temp_stim = cond.tinmeg3;
    end

    run('clean.m')

end

%% Run ICA on downsampled data to remove heartbeat artefacts.

%Run preprocessing for each subject
for i = 1:numel(sub_date.ID)

    inpath = ['../processed_data/preprocessed/' 'ID' sub_date.ID{i} '/'];
    outdir = ['../processed_data/ICA/' 'ID' sub_date.ID{i} '/'];
    figoutdir = ['../processed_data/ICA/figures/ID' sub_date.ID{i} '/'];

    %Check if dir for sub exist and skip
    if exist(outdir, 'file');
        warning(['ICA already exist for subject: ' sub_date.ID{i}])
    continue
    end

    if sub_date.Exp{i} == 'tinmeg1'
       temp_cond = cond.tinmeg1.all;
       temp_stim = cond.tinmeg1;
    elseif sub_date.Exp{i} == 'tinmeg2'
       temp_cond = cond.tinmeg2.all;
       temp_stim = cond.tinmeg2;
    elseif sub_date.Exp{i} == 'tinmeg3'
       temp_cond = cond.tinmeg3.all;
       temp_stim = cond.tinmeg3;
    end

    run('ICA.m')

end

%% Process EOG

for i = 1:numel(sub_date.ID)

    outdir = ['../processed_data/timelockeds/' 'ID' sub_date.ID{i} '/EOG/'];
    figoutdir = ['../processed_data/timelockeds/figures/' 'ID' sub_date.ID{i} '/'];

    %Check if dir for sub exist and skip or create
    if exist(outdir, 'file');
        warning(['EOG folder already exist for subject: ' sub_date.ID{i}])
    continue
    elseif ~exist(outdir, 'file');
        mkdir(outdir);
    end
    
    %What experiement version
    exp_ver = sub_date.Exp{i};

    if sub_date.Exp{i} == 'tinmeg1'
       temp_cond = cond.tinmeg1.all;
       temp_stim = cond.tinmeg1;
    elseif sub_date.Exp{i} == 'tinmeg2'
       temp_cond = cond.tinmeg2.all;
       temp_stim = cond.tinmeg2;
    elseif sub_date.Exp{i} == 'tinmeg3'
       temp_cond = cond.tinmeg3.all;
       temp_stim = cond.tinmeg3;
    end

    run('process_EOG.m')

end

%% MRI processing - Step 1: Require some manual input of fiducials and confirming coordsys

for i = 1:numel(sub_date.ID);
    
    mri_path_in = [mri_data_path 'NatMEG_' sub_date.MRI_ID{i}];
    mri_path_out = ['../processed_data/MRI/ID' sub_date.ID{i} '/'];

    mri_seq = sub_date.MRI_seq{i}; %What MRI sequence to read

    %Continue if MRI_ID is empty for subject
    if (sub_date.MRI_ID{i} == "");
        warning(['NO MRI_ID exist for subject: ID' sub_date.ID{i}])
        continue
    end
    
    %Check if realigned mri exist for subject or create folder
    if exist([mri_path_out 'mri_realigned.mat'], 'file')
        warning(['Output mri_realigned exist for subject: ID' sub_date.ID{i}])
        continue
    elseif ~exist(mri_path_out, 'file')
        mkdir(mri_path_out);
    end
    
    run('MRI_processing1.m');
    
end

clear('mri_coordsys', 'subpath', 'mri_seq', 'outdir', 'mri_realigned_1', 'mri_realigned_vol_2', 'mri_realigned_vol_3', 'submri', 'sensshape', 'headshape', 'dcmfile', 'dicom_path', 'MEGfile', 'mri_path_in', 'mri_path_out', 'i');

%% MRI processing - Step 2: Create headmodel from MRI

for i = 1:numel(sub_date.ID);
    
    mri_path_in = ['../processed_data/MRI/ID' sub_date.ID{i} '/'];
    mri_path_out = ['../processed_data/MRI/ID' sub_date.ID{i} '/'];

    %Check if headmodel for subject or create folder
    if exist([mri_path_out 'headmodel.mat'], 'file')
        warning(['Headmodel exist for subject: ID' sub_date.ID{i}])
        continue
    elseif ~exist([mri_path_in 'mri_realigned.mat'], 'file')
        warning(['No mri_realigned found for subject: ID' sub_date.ID{i}])
        continue
    elseif ~exist(mri_path_out, 'file')
        mkdir(mri_path_out);
    end
    
    run('MRI_processing2.m');
    
end

%% Timelockeds

for i = 1:numel(sub_date.ID)

    inpath = ['../processed_data/ICA/' 'ID' sub_date.ID{i} '/'];
    outdir = ['../processed_data/timelockeds/' 'ID' sub_date.ID{i} '/'];
    figoutdir = ['../processed_data/timelockeds/figures/' 'ID' sub_date.ID{i} '/'];

    %What experiement version
    exp_ver = sub_date.Exp{i};

    if sub_date.Exp{i} == 'tinmeg1'
       temp_cond = cond.tinmeg1.all;
       temp_stim = cond.tinmeg1;
    elseif sub_date.Exp{i} == 'tinmeg2'
       temp_cond = cond.tinmeg2.all;
       temp_stim = cond.tinmeg2;
    elseif sub_date.Exp{i} == 'tinmeg3'
       temp_cond = cond.tinmeg3.all;
       temp_stim = cond.tinmeg3;
    end

    run('process_tlks.m')

end

%% Gather timelockeds

%Load if not already in workspace
if exist('tlk_all_sub','var') == 0;
    tlk_all_sub = load('../processed_data/timelockeds/aggregated/tlk_all_sub_cmb.mat');
    tlk_all_sub = tlk_all_sub.tlk_all_sub;
end

%create empty struct if file doesnd exist

%For each subject
for i = 1:numel(sub_date.ID)

    %Get experiment version
    temp_exp = (sub_date.Exp{i});
    
    IDs = find(~cellfun(@isempty, tlk_all_sub.(temp_exp).ID));
    
    %Find first empty cell for ID
    empty_IDx = 1 + IDs(end);
    
    %Check if subject is already in struct and skip
    if any(ismember(tlk_all_sub.(temp_exp).ID, sub_date.ID{i}))
        warning([sub_date.ID{i} ' is already in tlk structure'])
        continue
    else
        run('gather_tlks.m')
    
    %if subject already in struct
    end

%For subject
end

save('../processed_data/timelockeds/aggregated/tlk_all_sub_cmb.mat', 'tlk_all_sub');

%% Gather EOG

%Load if not already in workspace
if exist('EOG_all_sub','var') == 0;

    %Create empty if file doesnt exist
    if exist('../processed_data/timelockeds/aggregated/EOG_all_sub.mat', 'file') == 0;
        
        %Create empty structure
        EOG_all_sub = struct();
    
    %else load
    else
        EOG_all_sub = load('../processed_data/timelockeds/aggregated/EOG_all_sub.mat');
        EOG_all_sub = EOG_all_sub.EOG_all_sub;
    end
end

%For each subject
for i = 1:numel(sub_date.ID)

    %Get experiment version
    temp_exp = (sub_date.Exp{i});

    %If no field for experiment, create it
    if ~any(ismember(fieldnames(EOG_all_sub), temp_exp))
       EOG_all_sub.(temp_exp).ID = {};
    end
    
    %if ID array is empty, write to first line
    if isempty(EOG_all_sub.(temp_exp).ID)
        empty_IDx = 1;

    %else, find the first empty cell
    else
        IDs = find(~cellfun(@isempty, EOG_all_sub.(temp_exp).ID));
        empty_IDx = 1 + IDs(end);
    end
    
    %Check if subject is already in struct and skip
    if any(ismember(EOG_all_sub.(temp_exp).ID, sub_date.ID{i}))
        warning([sub_date.ID{i} ' is already in tlk structure'])
        continue
    else
        run('gather_EOG.m')
    
    %if subject already in struct
    end

%For subject
end

save('../processed_data/timelockeds/aggregated/EOG_all_sub.mat', 'EOG_all_sub');

%% Clean EOG

%Load EOG data to clean up
if exist('EOG_all_sub','var') == 0;
    EOG_all_sub = load('../processed_data/timelockeds/aggregated/EOG_all_sub.mat');
    EOG_all_sub = EOG_all_sub.EOG_all_sub;
end

%Load if not already in workspace
if exist('EOG_clean','var') == 0;

    %Create empty if file doesnt exist
    if exist('../processed_data/timelockeds/aggregated/EOG_clean.mat', 'file') == 0;
        
        %Create empty structure
        EOG_clean = struct();
    
    %else load
    else
        EOG_clean = load('../processed_data/timelockeds/aggregated/EOG_clean.mat');
        EOG_clean = EOG_clean.EOG_clean;
    end
end

%For each subject
for i = 1%:numel(sub_date.ID)

    %Get experiment version
    temp_exp = (sub_date.Exp{i});

    %If no field for experiment, create it
    if ~any(ismember(fieldnames(EOG_clean), temp_exp))
       EOG_clean.(temp_exp).ID = {};
    end
    
    %if ID array is empty, write to first line
    if isempty(EOG_clean.(temp_exp).ID)
        empty_IDx = 1;

    %else, find the first empty cell
    else
        IDs = find(~cellfun(@isempty, EOG_clean.(temp_exp).ID));
        empty_IDx = 1 + IDs(end);
    end
    
    %Check if subject is already in struct and skip
    if any(ismember(EOG_clean.(temp_exp).ID, sub_date.ID{i}))
        warning([sub_date.ID{i} ' is already in EOG_clean structure'])
        continue
    else
        run('clean_EOG.m')
    
    %if subject already in struct
    end

%For subject
end

save('../processed_data/timelockeds/aggregated/EOG_clean.mat', 'EOG_clean');
