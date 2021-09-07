%% Read MRI

%https://github.com/natmegsweden/meeg_course/blob/master/tutorial_03_prepare_mri.md

%% Read and process DICOM/nifti MRI for subjects with no headmodel

%To do:
%Add title of subject ID in figure for checks

%% Check for dicoms, read, coordinate, realign and save for all participants

%i = sub_date.ID(1);

%Check for DICOM/nifti and read_mri
if exist(fullfile(mri_path_in, '/DICOM'))

dicom_path = [mri_data_path 'NatMEG_' sub_date.ID{i} '/DICOM/'];
dcmfile = dir(dicom_path);
dcmfile = [dcmfile(3).name]; %Grab 3rd to avoid folders

submri = ft_read_mri([dicom_path dcmfile]);

elseif exist(fullfile(mri_path_in, '/nifti'))
submri = ft_read_mri([mri_path_in '/nifti/NatMEG_' sub_date.ID{i} '_T1w.nii.gz']);

end
%Inspection plot
%ft_sourceplot([], submri)

%Determine coordsys directions
warning(['Now displaying coordsys for subject ID ' sub_date.ID{i}])
mri_coordsys = ft_determine_coordsys(submri);

save([mri_path_out 'mri_coordsys'], 'mri_coordsys');

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';

%Provide n, r, l and positive Z (e.g. top of scalp), q to save and quit
warning(['Now displaying volume realign for subject ID ' sub_date.ID{i}])
mri_realigned_1 = ft_volumerealign(cfg, mri_coordsys);

%save([mri_path_out 'mri_realigned_1'], 'mri_realigned_1');

%Read MEG data for subject, takes first available fif
subpath = [meg_data_path 'NatMEG_' sub_date.ID{i} '/' sub_date.date{i} '/'];
fnames = find_files(subpath, {'tinmeg1', 'tsss'}, 'ds');
fnames = fnames{1};

MEGfile = fullfile(subpath, fnames);

headshape = ft_read_headshape(MEGfile);
sensshape = ft_read_sens(MEGfile);

save([mri_path_out 'headshape'], 'headshape');
save([mri_path_out 'sensshape'], 'sensshape');

%plot headshape points
% figure
% ft_plot_headshape(headshape);
% ft_plot_sens(sensshape);

%Check volume_realign is correct (nasion)
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

save([mri_path_out 'mri_realigned'], 'mri_realigned_vol_3');
