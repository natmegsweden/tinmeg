%% Read MRI

%https://github.com/natmegsweden/meeg_course/blob/master/tutorial_03_prepare_mri.md

%% Filepaths
meg_data_path = '/archive/20061_tinnitus/MEG/';

mri_data_path = '../MRI/';

%load subjects list
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

%% Read and process DICOM/nifti MRI for subjects with no headmodel

%To do:
%Add title of subject ID in figure for checks

%Check for dicoms, read, coordinate, realign and save for all participants
for i = 1:2%:length(sub_date.ID);
    
    sub_mri_path = ['../MRI/' 'NatMEG_' char(sub_date{i,1})];
    fname = ['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_MEG_headmodel.mat'];
    
    %Check if headmodel exist for subject
    if exist(fname, 'file')
    warning(['Output' fname ' exist for subject ' char(sub_date{i,1})])
    continue
    end
    
    %Check for DICOM/nifti and read_mri
    if exist(fullfile(sub_mri_path, '/DICOM'))
    
    dicom_path = [mri_data_path 'NatMEG_' sub_date.ID{i} '/DICOM/'];
    dcmfile = dir(dicom_path);
    dcmfile = [dcmfile(3).name]; %Grab 3rd to avoid folders
    
    submri = ft_read_mri([dicom_path dcmfile])
    
    elseif exist(fullfile(sub_mri_path, '/nifti'))
    submri = ft_read_mri([sub_mri_path '/nifti/NatMEG_' char(sub_date{i,1}) '_T1w.nii.gz'])
    
    end
    %Inspection plot
    %ft_sourceplot([], submri)
    
    %Determine coordsys directions
    mri_coordsys = ft_determine_coordsys(submri);
    
    %save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_coordsys'], 'mri_coordsys');
    
    cfg = [];
    cfg.method = 'interactive';
    cfg.coordsys = 'neuromag';
    
    %Provide n, r, l and positive Z (e.g. top of scalp), q to save and quit
    mri_realigned_1 = ft_volumerealign(cfg, mri_coordsys);
    
    %save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_realigned_1'], 'mri_realigned_1');
    
    %Read MEG data for subject, takes first available fif
    subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
    fnames = find_files(subpath, {'tinmeg1', 'tsss'}, 'ds');
    fnames = fnames{1};
    
    MEGfile = fullfile(subpath, fnames);

    headshape = ft_read_headshape(MEGfile);
    sensshape = ft_read_sens(MEGfile);
    
%     %plot headshape points
%     figure
%     ft_plot_headshape(headshape);
%     ft_plot_sens(sensshape);
       
    %Check volume_realign is correct (nasion) if not
    cfg = [];
    cfg.method              = 'headshape';
    cfg.headshape.headshape = headshape;
    cfg.headshape.icp       = 'yes';
    cfg.coordsys            = 'neuromag';

    mri_realigned_vol_2 = ft_volumerealign(cfg, mri_realigned_1);
    
    %Check volume_realign is correct after processing
    cfg = [];
    cfg.method              = 'headshape';
    cfg.headshape.headshape = headshape;
    cfg.headshape.icp       = 'no'; %Do not fit again.
    cfg.coordsys            = 'neuromag';

    mri_realigned_vol_3 = ft_volumerealign(cfg, mri_realigned_vol_2);
    
    %save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_realigned_3'], 'mri_realigned_vol_3');
       
    %Reslice MRI
    cfg = [];
    cfg.resolution = 1;

    mri_resliced = ft_volumereslice(cfg, mri_realigned_vol_3);
    
    %Convert units
    mri_resliced = ft_convert_units(mri_resliced, 'cm');
    
    %save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_resliced'], 'mri_resliced');
    
    %Segment MRI
    cfg = [];
    cfg.output = {'brain' 'skull' 'scalp'};
    
    mri_segmented = ft_volumesegment(cfg, mri_resliced);

    %save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_segmented'], 'mri_segmented');
    
    %Correct compartments
    binary_brain = mri_segmented.brain;
    binary_skull = mri_segmented.skull | binary_brain;
    binary_scalp = mri_segmented.scalp | binary_brain | binary_skull;

    % use boolean logic together with IMERODE
    binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
    binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull
    
    % Copy MRI
    mri_segmented_2 = mri_segmented;                    % Copy stucture
    mri_segmented_2.anatomy = mri_resliced.anatomy;  % Copy anatomical data
  
    % insert the updated binary volumes, taking out the center part for skull and scalp
    mri_segmented_2.brain    = binary_brain;
    mri_segmented_2.skull    = binary_skull & ~binary_brain;
    mri_segmented_2.scalp    = binary_scalp & ~binary_brain & ~binary_skull;
    
    
%     % Plot segmentations
%     cfg = [];
%     cfg.funparameter = 'brain';
%     ft_sourceplot(cfg, mri_segmented_2);
% 
%     cfg.funparameter = 'skull';
%     ft_sourceplot(cfg, mri_segmented_2);
% 
%     cfg.funparameter = 'scalp';
%     ft_sourceplot(cfg, mri_segmented_2);

    %save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_segmented_2'], 'mri_segmented_2');


    %Create mesh for brain
    cfg = [];
    cfg.method = 'projectmesh';
    cfg.tissue = 'brain';
    cfg.numvertices = 3000;

    mesh_brain = ft_prepare_mesh(cfg, mri_segmented_2);
    
    %ft_plot_mesh(mesh_brain, 'facealpha', .5, 'edgealpha', 0.1)
    
    %save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mesh_brain'], 'mesh_brain');
    
    
    %Final head_model (only need brain for MEG);
    cfg = [];
    cfg.method = 'singleshell';

    headmodel_meg = ft_prepare_headmodel(cfg, mesh_brain);
    
    save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_MEG_headmodel'], 'headmodel_meg');
    
    %Final plot - aligned MEG
    figure; hold on
    ft_plot_sens(sensshape)
    ft_plot_headshape(headshape)
    ft_plot_headmodel(headmodel_meg)
    ft_plot_axes([], 'unit', 'cm');
    
%     May use this to pause loop while figure open, i.e to inspect plots
%     h = axes;
% 
%     while ishandle(h);
%     data = rand(1000,10);
%     plot(data); hold on
%     drawnow;    
%     end
    
end
