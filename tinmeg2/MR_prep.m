%% MR_prep

%% Check for dicoms, read, coordinate, realign and save for all participants

addpath('../')

mri_path_in = ['../../MRI/' 'NatMEG_' sub_date.MRI_ID{i}];
mri_path_out = ['../../mat_data/MRI_mat_tinmeg2/ID' sub_date.ID{i} '/'];

%Check for DICOM/nifti and read_mri
if exist(fullfile([mri_path_in, '/DICOM']))

    dicom_path = [mri_path_in '/DICOM/'];
    dcmfile = dir(dicom_path);
    dcmfile = [dcmfile(3).name]; %Grab 3rd to avoid folders
    submri = ft_read_mri([dicom_path dcmfile]);

elseif exist(fullfile([mri_path_in, '/nifti']))

    submri = ft_read_mri([mri_path_in '/nifti/NatMEG_' sub_date.MRI_ID{i} '_T1w.nii.gz']);

elseif exist(fullfile(mri_path_in, '/T1.mgz'))
    submri = ft_read_mri([mri_path_in '/T1.mgz']);
end

%Inspection plot
%ft_sourceplot([], submri)

%Determine coordsys directions
warning(['Now displaying coordsys for subject ID ' sub_date.MRI_ID{i}])
mri_coordsys = ft_determine_coordsys(submri);

close;

save([mri_path_out 'mri_coordsys.mat'], 'mri_coordsys');

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';

%Provide n, r, l and positive Z (e.g. top of scalp), q to save and quit
warning(['Now displaying volume realign for subject ID ' sub_date.ID{i}])
mri_realigned_1 = ft_volumerealign(cfg, mri_coordsys);

%save([mri_path_out 'mri_realigned_1'], 'mri_realigned_1');

%Read MEG data for subject, takes first available fif
subpath = [meg_data_path 'NatMEG_' sub_date.ID{i} '/' sub_date.date{i} '/'];
fnames = find_files(subpath, {'tinmeg2', 'tsss'}, 'ds');
fnames = fnames{1};

MEGfile = fullfile(subpath, fnames);

headshape = ft_read_headshape(MEGfile);
sensshape = ft_read_sens(MEGfile);
cd
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