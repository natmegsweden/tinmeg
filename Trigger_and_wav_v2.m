
%% Triggers structure and audiofile names

triggers = struct();

triggers.tin0_bkg0 = [0 1 2 4 8 16 32];
triggers.tin0_bkg3 = triggers.tin0_bkg0([1 4:7]) + 48;
triggers.tin0_bkg8 = triggers.tin0_bkg0([1 4:7]) + 32;

triggers.tin3_bkg0 = triggers.tin0_bkg0([1 4:7]) + 192;
triggers.tin3_bkg3 = triggers.tin0_bkg0([1 4:7]) + 240;
triggers.tin3_bkg8 = triggers.tin0_bkg0([1 4:7]) + 224;

triggers.tin8_bkg0 = triggers.tin0_bkg0([1 4:7]) + 128;
triggers.tin8_bkg3 = triggers.tin0_bkg0([1 4:7]) + 176;
triggers.tin8_bkg8 = triggers.tin0_bkg0([1 4:7]) + 160;

audiofiles = {'tin0_bkg0.wav',
              'tin0_bkg3.wav',
              'tin0_bkg8.wav',
              'tin3_bkg0.wav',
              'tin3_bkg3.wav',
              'tin3_bkg8.wav',
              'tin8_bkg0.wav',
              'tin8_bkg3.wav',
              'tin8_bkg8.wav'}

addpath('../../audio/');

conditions = fieldnames(triggers);

%% Loop over conditions in fif and extract triggers

%read in fif
%infile = '/home/nikedv/TinMEG1/headpos_output/NatMEG_0905/221006/nomeg_trigtest.fif'

% 5kHz recording alternative:
%infile = '/archive/20061_tinnitus/MEG/NatMEG_0905/221006/nomeg_tintest_ny.fif'

infile = '/archive/20061_tinnitus/MEG/NatMEG_0905/221010/nomeg_tintest.fif'

%cfg for read header and event
cfg                     = [];
cfg.dataset             = infile;
cfg.trialdef.prestim    = 1;
cfg.trialdef.poststim   = 1;
cfg.trialdef.eventtype = 'STI101';

% read the header information and the events from the data
hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset, 'detectflank', 'up');

%Pivot events and convert to cell
event_cell = struct2cell(event)';

%Only interested in subset of channel STI101
event_101 = strcmp(event_cell, 'STI101');
event_cell = event_cell(event_101, :);

%Onlys interested in Trigger values and Onsets - convert to matrix.
trigsamples = cell2mat([event_cell(:,3) event_cell(:,2)]);

for i = 1:(numel(conditions))
    
    %Read in audio file
%     audf = audioread(['../../audio/' audiofiles{i}]);

    %Convert to "mono"
%     audf = audf(:,1);
    
    %Get only triggers relevant for condition
    condtrigs = trigsamples(ismember(trigsamples(:,1), triggers.(conditions{i})),:);

    %Normalize to start at time zero (i.e. trigger for soundfile start)
    condtrigs(:,2) = condtrigs(:,2) - condtrigs(1,2);

    %convert samples from 1k to 44.1kHz fs and round to integer
    condtrigs(:,2) = round(condtrigs(:,2)*(44100/hdr.Fs),0);

    %Convert to milliseconds and write the first 21 triggers (one trial block) to cell
    %Note first condition "tin0_bkg0" has filestart trigger = 0, ignore for evaluation
    for ii = 1:11

    fiftrig_cum(ii,i) = round(condtrigs(ii,2)/44.100, 0);

    end
    
end

varnames = conditions';
%writetable(array2table(fiftrig_cum, 'VariableNames', varnames),'fiftriggers2.txt', 'Delimiter', '\t')
