
load(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_realigned']);

%Reslice MRI
cfg = [];
cfg.resolution = 1;

mri_resliced = ft_volumereslice(cfg, mri_realigned_vol_3);

%Convert units
mri_resliced = ft_convert_units(mri_resliced, 'cm');

save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_resliced'], 'mri_resliced');

%Segment MRI
cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};

mri_segmented = ft_volumesegment(cfg, mri_resliced);

save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_segmented'], 'mri_segmented');

%Correct compartments
binary_brain = mri_segmented.brain;
binary_skull = mri_segmented.skull | binary_brain;
binary_scalp = mri_segmented.scalp | binary_brain | binary_skull;

% use boolean logic together with IMERODE
binary_skull = binary_skull & imerode(binary_scalp, strel_bol(2)); % fully contained inside eroded scalp
binary_brain = binary_brain & imerode(binary_skull, strel_bol(2)); % fully contained inside eroded skull

% Copy MRI
mri_segmented_2 = mri_segmented;                 % Copy stucture
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

save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mri_segmented_2'], 'mri_segmented_2');

%Create mesh for brain
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;

mesh_brain = ft_prepare_mesh(cfg, mri_segmented_2);

%ft_plot_mesh(mesh_brain, 'facealpha', .5, 'edgealpha', 0.1)

save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_mesh_brain'], 'mesh_brain');


%Final head_model (only need brain for MEG);
cfg = [];
cfg.method = 'singleshell';

headmodel_meg = ft_prepare_headmodel(cfg, mesh_brain);

save(['../mat_data/MRI_mat/' 'ID' char(sub_date{i,1}) '_MEG_headmodel'], 'headmodel_meg');

