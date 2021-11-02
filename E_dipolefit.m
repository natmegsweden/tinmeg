
%https://github.com/natmegsweden/meeg_course/blob/master/tutorial_04a_dipole_fitting.md

for i = 1%:length(sub_date.ID)

    %load leadfield (from E_beamformer) - Subject to change
    %leadfield = load(['../mat_data/source_reconstruction/' 'ID' sub_date.ID{i} '/leadfield.mat']);
    %leadfield = leadfield.leadfield;
    
    %load headmodel (from D_MR_prep2)
    headmodel = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/meg_headmodel.mat']);
    headmodel = headmodel.headmodel_meg;
    
    %load reslice MRI for plot
    mri = load(['../mat_data/MRI_mat/ID' sub_date.ID{i} '/mri_resliced.mat']);
    mri = mri.mri_resliced;
    
    for ii = 1%:length(conditions)
        
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    trig = eval(['cond.' char(conditions(ii)) 'trig']);

        for iii = 5%1:nstim
        
        timelockeds = load(['../mat_data/timelockeds/ID' sub_date.ID{i} '/' (cond.(([conditions{ii} 'label'])){iii}) '_tlks.mat']);
        timelockeds = timelockeds.timelockeds;
        
        %Make leadfields for magnetometers
        cfg.senstype        = 'meg';
        cfg.grad            = timelockeds.grad;
        cfg.headmodel       = headmodel;
        cfg.channel         = 'megmag';
        
        leadfield_mag = ft_prepare_leadfield(cfg, timelockeds);
        
        cfg = [];
        cfg.gridsearch      = 'yes';            % search the grid for an optimal starting point
        cfg.numdipoles      = 2;                % N dipoles in model
        cfg.symmetry        = 'x';              % Leave empty for single dipole fit
        cfg.sourcemodel     = leadfield_mag;    % supply the grid
        cfg.headmodel       = headmodel;        % supply the headmodel
        cfg.dipfit.metric   = 'rv';             % the metric to minimize
        cfg.model           = 'regional';       % Assume that the dipole has a fixed position
        cfg.senstype        = 'meg';            % sensor type
        cfg.channel         = 'megmag';         % which channels to use
        cfg.nonlinear       = 'yes';            % do a non-linear search
        cfg.latency         = [0.050 0.150];    % specify the TOI
        cfg.backproject     = 'yes';            % Predict values from model

        %ft_dipolefitting from fieldtrip-20200224, error with reshape (#450) in current version.
        dipole_mag_PN1 = ft_dipolefitting(cfg, timelockeds);
        
            
        %For subject
        end
        
    %For stim    
    end
    
%For condition    
end
