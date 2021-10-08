%% Loop over dowsampled data in ft_rejectvisual where no cleaned output

%Check if ID is in trial-log and determine row in log to write to
if find(strcmp(['ID' sub_date.ID{i}], cleancondlog)) > 0;
   logheight = find(strcmp(['ID' sub_date.ID{i}], cleancondlog));
else

%Find height of trial-log and +1 for new row
logheight = size(cleancondlog);
logheight = logheight(1) + 1;

%Write ID to new row, column 1
cleancondlog{logheight,1} = ['ID' sub_date.ID{i}];
end

%For each condition loaded
for ii = 1:length(conditions);
fname = [char(conditions(ii)) '_ds_clean' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

    %check if file exist
    if exist(fpath, 'file')
    warning([char(conditions(ii)) ' for subject: ID' sub_date.ID{i} ' exist'])
    continue
    end

%if not exist load downsampled (_ds) mat-file
fname = [char(conditions(ii)) '_ds' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

res4mat_ds = load(fpath);
res4mat_ds = res4mat_ds.res4mat_ds

%ft_rejectvisual for condition
cfg = [];
cfg.method = 'summary';
cfg.keepchannel = 'yes';
cfg.channel = 'MEGMAG';
cfg.layout = 'neuromag306eeg_1005_natmeg.lay';

cleaned4mat = ft_rejectvisual(cfg,res4mat_ds);

cfg.channel = 'MEGGRAD';

cleaned4mat = ft_rejectvisual(cfg, cleaned4mat);

%new filename and save
fname = [char(conditions(ii)) '_ds_clean' '.mat'];
fpath = ['../mat_data/preprocessing/' 'ID' sub_date.ID{i} '/' fname];

save(fpath, 'cleaned4mat');

    %write n of trials log        
    %How many stimtypes for cond and what trigger values
    nstim = length(eval(['cond.' char(conditions(ii)) 'trig']));
    trigs = eval(['cond.' char(conditions(ii)) 'trig']);

    %Write log for n of trials in condition after cleaning
    for iii = 1:length(eval(['cond.' char(conditions(ii)) 'trig']))
    if strcmp(conditions(ii), 'PO60')
    cleancondlog(logheight,iii+1) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'PO70')
    cleancondlog(logheight,iii+7) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP60')
    cleancondlog(logheight,iii+12) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GP70')
    cleancondlog(logheight,iii+16) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
    elseif strcmp(conditions(ii), 'GO')
    cleancondlog(logheight,iii+20) = num2cell(sum(cleaned4mat.trialinfo(:) == trigs(iii)));
    end
    end

%Write log to csv
writetable(cell2table(cleancondlog), '../Analysis Output/n_cond_clean.csv', 'WriteVariableNames', false) %Write log

%for conditions
end
