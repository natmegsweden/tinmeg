
%Check for DICOM folder/.nii.gz and read_mri
if exist([mri_path_in '/0000000' mri_seq '/'], 'dir')

dicom_path = [mri_path_in '/0000000' mri_seq '/'];
dcmfile = dir(dicom_path);
dcmfile = [dcmfile(3).name]; %Grab 3rd to avoid folders ft_read_mri finds all in series

submri = ft_read_mri([dicom_path dcmfile]);

%If no DICOM check for nii.gz
elseif exist([mri_path_in '/NatMEG_' sub_date.ID{i} '_T1w.nii.gz'])
submri = ft_read_mri([mri_path_in '/NatMEG_' sub_date.ID{i} '_T1w.nii.gz']);

end

%Inspection plot
%ft_sourceplot([], submri)

%Determine coordsys directions
warning(['Now displaying coordsys for subject ID ' sub_date.ID{i}])
mri_coordsys = ft_determine_coordsys(submri);

save([mri_path_out 'mri_coordsys'], 'mri_coordsys'); clear submri;

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';

%Provide n, r, l and positive Z (e.g. top of scalp), q to save and quit
warning(['Now displaying volume realign for subject ID ' sub_date.ID{i}])
mri_realigned_1 = ft_volumerealign(cfg, mri_coordsys);

%save([mri_path_out 'mri_realigned_1'], 'mri_realigned_1');

%Read MEG data for subject, takes first available fif
if sub_date.Exp{i} == 'TinMEG1'
    subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
    fnames = find_files(subpath, {'tinmeg1', 'tsss'}, 'ds');
elseif sub_date.Exp{i} == 'TinMEG2'
    subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
    fnames = find_files(subpath, {'tinmeg2', 'tsss'}, 'ds');
elseif sub_date.Exp{i} == 'TinMEG3'
    subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/meg/'];
    fnames = find_files(subpath, {'tinmeg3', 'tsss'}, 'ds');
end

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