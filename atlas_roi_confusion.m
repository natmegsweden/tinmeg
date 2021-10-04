%% ROI from atlas?

%loads as "pow_diff" this data is power difference from source reconstruction.
load('../mat_data/stats/PO60_90_pow_diff.mat');

aal_atlas = ft_read_atlas('../../fieldtrip-20210311/template/atlas/aal/ROI_MNI_V4.nii', 'unit', 'mm')

cfg = [];
cfg.parameter = 'tissue';
cfg.interpmethod = 'nearest';
int_aal_tempgrid = ft_sourceinterpolate(cfg, aal_atlas, template_grid);

int_aal_tempgrid.tissuelabel = aal_atlas.tissuelabel;

atlas_grid = ft_checkdata(int_aal_tempgrid, 'datatype', 'source');
atlas_grid.inside = template_grid.inside;

pow_diff.tissue = atlas_grid.tissue;
pow_diff.tissuelabel = atlas_grid.tissuelabel;

%loads as "mri"
load standard_mri

cfg = [];
cfg.parameter = 'tissue';
cfg.interpmethod = 'nearest';
int_powdiff_mri = ft_sourceinterpolate(cfg, pow_diff, mri);

labs = [81, 82];
atlas_grid.tissuelabel(labs)

int_powdiff_mri.tissue(~ismember(int_powdiff_mri.tissue, labs)) == 0;
int_powdiff_mri.tissue(ismember(int_powdiff_mri.tissue, labs)) == 1;

int_powdiff_mri.tissuelabel = {'ROI'};

cfg = [];
cfg.funparameter = 'tissue';
cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, int_powdiff_mri);

%% Or interpolate pow and parcellate to atlas?

%loads as int_powdiff, this data is power difference interpolated to standard mri.
load('../mat_data/stats/PO60_90_int_powdiff.mat');

cfg = [];

parcel = ft_sourceparcellate(cfg, int_powdiff, aal_atlas);

