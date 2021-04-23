%% Read MRI

%% Filepaths
meg_data_path = '/archive/20061_tinnitus/MEG/';

mri_data_path = '../MRI/'

%load subjects list
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

%% Read dicoms and realign

%To do
%ADD IF EXIST
%Add display of subject ID in figure for checks

%Check for dicoms, read, coordinate, realign and save for all participants
for i = 2%:length(sub_date.ID);
    
    dicom_path = [mri_data_path 'NatMEG_' sub_date.ID{i} '/DICOM/']
    dcmfile = dir(dicom_path);
    dcmfile = [dcmfile(3).name] %Grab 3rd to avoid folders
    
    submri = ft_read_mri([dicom_path dcmfile])
    
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
    fnames = fnames{1}
    
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
    
    save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_realigned_3'], 'mri_realigned_vol_3');
       
%     May use this to pause loop while figure open, i.e to inspect plots
%     h = axes;
% 
%     while ishandle(h);
%     data = rand(1000,10);
%     plot(data); hold on
%     drawnow;    
%     end
    
end
