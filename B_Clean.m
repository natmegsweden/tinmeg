%% Load downsampled data and clean with ft_rejectvisual


%% Read in subjects and conditions

% Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');
disp(['Number of subjects in table is ' num2str(height(sub_date))])

% Read conditions, triggers and their labels, to clean
load('../mat_data/conditions.mat');
load('../mat_data/cond.mat');

%% Loop over ft_rejectvisual where no cleaned output

%create log for n of trials in cleaned data
%first row are labels, first column IDs
if exist('../Analysis Output/n_cond_clean.csv', 'file');
    cleancondlog = readtable('../Analysis Output/n_cond_clean.csv', 'ReadVariableNames', false);
    cleancondlog = table2cell(cleancondlog);
else 
    %NB 22 columns hardcoded
    cleancondlog = ['ID' cond.PO60label cond.PO70label cond.GP60label cond.GP70label cond.GOlabel];
end

%For each subject in subject table
for i = 1:length(sub_date.ID);

    %Check if ID is in trial-log and determine row in log to write to
    if find(strcmp(['ID' char(sub_date.ID(i))], cleancondlog)) > 0;
       logheight = find(strcmp(['ID' char(sub_date.ID(i))], cleancondlog));
    else
        
    %Find height of trial-log and +1 for new row
    logheight = size(cleancondlog);
    logheight = logheight(1) + 1;
  
    %Write ID to new row, column 1
    cleancondlog{logheight,1} = ['ID' char(sub_date.ID(i))];
    end
    
    %For each condition loaded
    for ii = 1:length(conditions);
    fname = ['ID' char(sub_date.ID(i)) '_' char(conditions(ii)) '_ds_clean' '.mat'];
    fpath = ['../mat_data/' fname];

        %check if file exist
        if exist(fpath, 'file')
        warning([char(conditions(ii)) ' for subject: ID' char(sub_date.ID(i)) 'exist'])
        continue
        end
    
    %if not exist load downsampled (_ds) mat-file
    fname = ['ID' char(sub_date.ID(i)) '_' char(conditions(ii)) '_ds' '.mat']
    fpath = ['../mat_data/' fname]

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
    fname = ['ID' char(sub_date.ID(i)) '_' char(conditions(ii)) '_ds_clean' '.mat'];
    fpath = ['../mat_data/' fname];
    
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
    
%for subjects
end
