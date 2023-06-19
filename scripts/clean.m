%% Use Fieldtrip ft_rejectvisual to remove high variance trials/channels

%For each condition
for ii = 1:numel(temp_cond);
    fname = [temp_cond{ii} '_cl' '.mat'];
    fpath = ['../processed_data/preprocessed/' 'ID' sub_date.ID{i} '/' fname];
    
    %check if file exist
    if exist(fpath, 'file')
    warning([temp_cond{ii} ' for subject: ID' sub_date.ID{i} ' exist'])
    continue
    end
    
    %if not exist load downsampled (_ds) mat-file
    fname_in = [temp_cond{ii} '_ds' '.mat'];
    fpath_in = ['../processed_data/preprocessed/' 'ID' sub_date.ID{i} '/' fname_in];
    
    preproc_ds = load(fpath_in);
    preproc_ds = preproc_ds.preproc_ds;

    %ft_rejectvisual for condition
    cfg = [];
    cfg.method = 'summary';
    cfg.keepchannel = 'yes';
    cfg.channel = 'MEGMAG';
    cfg.layout = 'neuromag306all.lay';
    
    preproc_cl = ft_rejectvisual(cfg,preproc_ds);
    
    cfg.channel = 'MEGGRAD';
    
    preproc_cl = ft_rejectvisual(cfg, preproc_cl);
    
    %Save clean data
    save(fpath, 'preproc_cl');

    clear preproc_ds preproc_cl ii fname fname_in fpath fna_in fpath_in
    
%for conditions
end