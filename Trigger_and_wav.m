
%Confirming trigger timing to stimulation sound file for GP_ISI60

ISI60 = audioread('/archive/20061_tinnitus/Wav-filer/Audio Tracks/A_C60_i0.wav');

%Convert to "mono"
ISI60 = ISI60(:,1);

%Figure expect three triggers
%ISI60 start trigger        49729
%ISI60 gap-onset trigger    49730
%ISI60 pulse-onset trigger  49736

%ISI0 start                 49793
%ISI0 gap-onset trigger     49794
%ISI0 pulse-onset trigger   49800


%random example file
infile = '/archive/20061_tinnitus/MEG/NatMEG_0697/210208/tinmeg1-2_mc_avgtrans_tsss_corr98.fif'

%% from STI016fix to read triggers

cfg                     = [];   
cfg.dataset             = infile;
cfg.trialdef.prestim    = 1;
cfg.trialdef.poststim   = 1;
cfg.trialdef.eventvalue = [49793 49794 49800];
%cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';
%cfg                     = ft_definetrial(cfg);
%data                    = ft_preprocessing(cfg);

% Get only specific event type
cfg.trialdef.eventtype = ft_getopt(cfg.trialdef, 'eventtype', 'STI101');
cfg.trialdef.eventvalue = ft_getopt(cfg.trialdef, 'eventvalue', []);

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);

% Read trigger channels
chanindx = find(not(cellfun('isempty', strfind(hdr.label,'STI0'))));

event = ft_read_event(cfg.dataset, 'chanindx', chanindx);

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

%%

%Take samplepoints identified for triggers in trial function
trigsamples = trl(1:3,1:2);

%Normalize to start at time zero
trigsamples = trigsamples - trigsamples(1,1);

%convert samples from 5k to 44.1kHz fs and round to integer
trigsamples = round(trigsamples*(44100/5000),0);

%Figure of soundfile and triggers
figure
hold on
plot(ISI60)

plot([trigsamples(1,1) trigsamples(1,1)], [0 0.25], 'r', 'LineWidth', 2)
plot([trigsamples(2,1) trigsamples(2,1)], [0 0.25], 'r', 'LineWidth', 2)
plot([trigsamples(3,1) trigsamples(3,1)], [0 0.25], 'r', 'LineWidth', 2)

xlim([0 44100]);

hold off