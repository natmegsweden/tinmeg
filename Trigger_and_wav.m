
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

infile = '/archive/20061_tinnitus/MEG/NatMEG_0905/221006/nomeg_tintest_ny.fif'

cfg = [];

dat = ft_read_data(cfg, infile);

for i = 1:(numel(conditions))

    audf = audioread(['../../audio/' audiofiles{i}]);

    %Convert to "mono"
    audf = audf(:,1);
    
    event = ft_read_event(infile, 'detectflank', 'up');
    
    cfg =[];
    cfg.event = event;
    cfg.dataset = infile;
    cfg.trialfun = 'ft_trialfun_general';
    cfg.trialdef.eventtype = 'STI101';
    cfg.trialdef.prestim = 2.5;
    cfg.trialdef.poststim = 2.5;
    cfg.trialdef.eventvalue = triggers.(conditions{i});
    
    dat = ft_read_data(cfg.dataset);
    
    %from STI016fix to read triggers
    cfg                     = [];
    cfg.dataset             = infile;
    cfg.trialdef.prestim    = 1;
    cfg.trialdef.poststim   = 1;
    cfg.trialdef.eventvalue = triggers.(conditions{i});
    %cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';
    %cfg                     = ft_definetrial(cfg);
    %data                    = ft_preprocessing(cfg);

    % Get only specific event type
    cfg.trialdef.eventtype = ft_getopt(cfg.trialdef, 'eventtype', 'STI101');
    cfg.trialdef.eventvalue = ft_getopt(cfg.trialdef, 'eventvalue', []);

    % read the header information and the events from the data
    hdr = ft_read_header(cfg.dataset);

    % Read trigger channels
    chanindx = find(not(cellfun('isempty', strfind(hdr.label,'STI0'))));

    %event = ft_read_event(cfg.dataset, 'chanindx', chanindx, 'detectflank', 'up');
    
    event = ft_read_event(cfg.dataset, 'detectflank', 'up');
    
    event_cell = struct2cell(event);

    % Manually make combined trigger channel
    dat = ft_read_data(cfg.dataset , 'chanindx', chanindx); % Read the trigger data
    allsti = dat(1:16,:)==5;
    trig = allsti(1,:)*2^0 + ...
           allsti(2,:)*2^1 + ...
           allsti(3,:)*2^2 + ...
           allsti(4,:)*2^3 + ...
           allsti(5,:)*2^4 + ...
           allsti(6,:)*2^5 + ...
           allsti(7,:)*2^6 + ...
           allsti(8,:)*2^7 + ...
           allsti(9,:)*2^8 + ...
           allsti(10,:)*2^9 + ...
           allsti(11,:)*2^10 + ...
           allsti(12,:)*2^11 + ...
           allsti(13,:)*2^12 + ...
           allsti(14,:)*2^13 + ...        
           allsti(15,:)*2^14 + ...        
           allsti(16,:)*2^15;

    % event = struct();
    for j=find(diff([0 trig])>0)
        event(end+1).type   = 'STI101';
        event(end  ).sample = j;          % assign the sample at which the trigger has gone up
        event(end  ).value  = trig(j);    % assign the trigger value just _after_ going up
    end

    %
    value  = [event.value]';
    sample = [event.sample]';

    %
    pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
    posttrig =  round(cfg.trialdef.poststim * hdr.Fs);

    if isempty(cfg.trialdef.eventvalue)
        cfg.trialdef.eventvalue = unique(value);
    end

    % look for the triggers
    trl = [];
    for j = 1:length(value)
        if strcmp(cfg.trialdef.eventtype, event(j).type)
            trg = value(j);
            if any(cfg.trialdef.eventvalue == trg)
                trlbegin = sample(j) + pretrig;       
                trlend   = sample(j) + posttrig;       
                offset   = pretrig;
                newtrl   = [trlbegin trlend offset, trg];
                trl      = [trl; newtrl];
            end
        end
    end

    %Take samplepoints identified for triggers in trial function
    trigsamples = [trl(:,1) trl(:,2)];

    %Normalize to start at time zero
    trigsamples = trigsamples - trigsamples(1,1);

    %convert samples from 1k to 44.1kHz fs and round to integer
    trigsamples = round(trigsamples*(44100/hdr.Fs),0);

    %Convert to milliseconds and write the first 21 triggers to cell
    for ii = 1:21

    fiftrig_cum(ii,i) = round(trigsamples(ii,1)/44.100, 0);

    end
    
end
varnames = conditions';
%writetable(array2table(fiftrig_cum, 'VariableNames', varnames),'fiftriggers2.txt', 'Delimiter', '\t')

%%

%Figure of soundfile and triggers
figure
hold on
plot(audf)
for i = 1:21
    plot([trigsamples(i,1) trigsamples(i,1)], [0 0.25], 'r', 'LineWidth', 2)
end
hold off
