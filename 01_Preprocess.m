%% Read and arrange fif-files for TinMEG1
%
%
%

%% Read subjects and dates

% Specify MEG data path
meg_data_path = '/archive/20061_tinnitus/MEG/';

%Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

disp(['Number of subjects in table is ' num2str(height(sub_date))])

%% Find relevant files for subject and create cell-array of file paths

% Create cell array for subjects files
fpaths = cell(1);

for i = 1:height(sub_date);

% Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
fnames = find_files(subpath, {'tinmeg1', 'tsss'}, 'ds');

fpaths{i,1} = char(sub_date.ID{i}); %Include ID for tracking

    for fileindex = 1:length(fnames);
        fpaths{i,1+fileindex} = [subpath char(fnames(fileindex))]; % NB! n of files differ between rows, some subjects have empty columns
    end
clear fnames subpath
end

writetable(cell2table(fpaths), '../Analysis Output/included_filepaths.csv') %Write log

%% Specify event triggers

%Trigger value (STI101) at pulse onset
eventsPO60 = [40968 36872 34824 33800 33288 33032];
eventsPO70 = [36876 34828 33804 33292 33036];
eventsGP60 = [49800 49736 49704 49688];
eventsGP70 = [49804 49740 49708 49692];

%Trigger value (STI101) at gap onset
eventsGO   = [16386 16390];

% Pulse trigger:
% B_C60_P70: 40968
% B_C60_P75: 36872
% B_C60_P80: 34824
% B_C60_P85: 33800
% B_C60_P90: 33288
% B_C60_P95: 33032

% B_C70_P75: 36876
% B_C70_P80: 34828
% B_C70_P85: 33804
% B_C70_P90: 33292
% B_C70_P95: 33036

%ISI pulse trigger (gaptrig = pulsetrig - 6)
% A_C60_i0: 49800
% A_C60_i60: 49736
% A_C60_i120: 49704
% A_C60_i240: 49688

% A_70_i0: 49804
% A_70_i60: 49740
% A_70_i120: 49708
% A_70_i240: 49692

%GO gap-trig:
% A_GO_60: 16386
% A_GO_70: 16390

%% Define trials and preprocess PO60

for i = 1:height(sub_date)

%Create filename
fname = ['ID' char(fpaths(i,1)) '_60PO_raw' '.mat']
fpath = ['../mat_data/' fname]

%check if file exist
if exist(fpath, 'file')
    warning(['Output exist for subject ' char(fpaths(i,1))])
    continue
end
    
% Define trials
cfg = [];

% NB! cellfun for cfg.dataset defines 2:highest populated column
cfg.dataset             = fpaths(i,2:max(find(~cellfun(@isempty,fpaths(i,:)))));
cfg.trialdef.prestim    = 0.35;        % seconds before trigger
cfg.trialdef.poststim   = 0.30;        % seconds after trigger
cfg.trialdef.eventvalue = eventsPO60;
cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';

cfg = ft_definetrial(cfg);

%Remove trials from cfg.trl that have negative sample index for trial start
cfg.trl = cfg.trl(cfg.trl(:,1) >= 0,:);
  
%Remove trials from cfg.trl that have higher sample index than exist in file
cfg.trl = cfg.trl(cfg.trl(:,2) < max([cfg.event.sample]),:);

% preprocessing
cfg.demean     = 'no';
cfg.lpfilter   = 'no';
cfg.hpfilter   = 'no';
cfg.dftfilter  = 'no';
cfg.channel    = {'MEG', 'ECG', 'EOG'};

res4mat = ft_preprocessing(cfg);

%log n of trials per subject
%res4mat.trialinfo

save(fpath, 'res4mat')

%downsample and save
cfg = [];
cfg.resamplefs = 200;

res4mat_ds = ft_resampledata(cfg, res4mat);

fname = ['ID' char(fpaths(i,1)) '_60PO_ds' '.mat'];
fpath = ['../mat_data/' fname];

save(fpath, 'res4mat_ds')

%clear temp variables
clear res4mat res4mat_ds

end

%%
%n trials - wip
for i = 1:n_subjects
    tottrig(i) = sum(res.(['ID_' fpaths{i,1}]).PO60.trialinfo == eventsPO60(1))
end

%% Create structure to read mat-files to

res = struct();

vararray = {'PO60', 'PO70', 'GP60', 'GP70', 'GO'};

for i = 1:n_subjects
    for ii = 1:length(vararray)
        res.(['ID_' fpaths{i,1}]).(vararray{ii}) = 0 + ii; %Concatenate 'ID_' as int in struct not allowed
    end
end

clear('i', 'ii', 'vararray')

%res.(['ID_' fpaths{i,1}]).PO60

% Tutorial snippets:
%
% for i = 1:n_subjects
%     epochs = struct('PO60', 1, 'PO70', 2, 'GP60', 3, 'GP70', 4, 'GO', 5)
% end
% 
% a = {'see','why'};
% KPI = {'L','L2','L3'};
% S.(a{1}).(KPI{1}) = 5;
% S.see.L
