
%% Line out from AudioFile to MISC001

infile = '/archive/20061_tinnitus/MEG/NatMEG_0905/220328/bench_audiomisc.fif';

misc01 = ft_read_data(infile, 'chanindx', 1);

% read the header information and the events from the data
hdr = ft_read_header(infile);
event = ft_read_event(infile, 'detectflank', 'up');

plot(misc01);


%% Check Audiofile on/off trigger

infile = '/archive/20061_tinnitus/MEG/NatMEG_0905/220323/soundfile_STI009.fif'; %Trigger every 1000ms

%cfg for read header and event
cfg                     = [];
cfg.dataset             = infile;
cfg.trialdef.prestim    = 1;
cfg.trialdef.poststim   = 1;
cfg.trialdef.eventtype = 'STI101';

% read the header information and the events from the data
hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset, 'detectflank', 'both');

%Pivot events and convert to cell
event_cell = struct2cell(event)';

%Only interested in subset of channel STI101
event_101 = strcmp(event_cell, 'STI009_up');
event_cell_up = event_cell(event_101, :);

event_101 = strcmp(event_cell, 'STI009_down');
event_cell_down = event_cell(event_101, :);

%% Read benchmarkfiles

%read in fif
onesec = '/archive/20061_tinnitus/MEG/NatMEG_0905/220323/benchmark.fif'; %Trigger every 1000ms
tensec = '/archive/20061_tinnitus/MEG/NatMEG_0905/220323/benchmark_10s.fif'; %Trigger every 10 000ms
newpres = '/archive/20061_tinnitus/MEG/NatMEG_0905/220323/benchmark_newpres.fif'; %Trigger every 1000ms, new presentation version

files = {onesec, tensec, newpres};
filenicks = {'onesec', 'tensec', 'newpres'};

trigsamples = struct();

for i = 1:numel(files);
    %cfg for read header and event
    cfg                     = [];
    cfg.dataset             = files{i};
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

    %Onlys interested in Trigger values and Onsets - convert to matrix
    tempdat = cell2mat([event_cell(:,3) event_cell(:,2)]);
    
    %New column normalize first sample to 0
    tempdat(:,3) = tempdat(:,2) - tempdat(1,2)
    
    %Hardcode in expected sample n
    if strcmp(filenicks{i}, 'onesec') || strcmp(filenicks{i}, 'newpres');
        tempdat(:,4) = (0:5000:(numel(tempdat(:,3))-1)*5000)';
    elseif strcmp(filenicks{i}, 'tensec')
       tempdat(:,4) = (0:50000:(numel(tempdat(:,3))-1)*50000)';   
    end
    
    tempdat(:,5) = tempdat(:,3) - tempdat(:,4);
    
    %New column translate to millisecond
    tempdat(:,6) = round(tempdat(:,3)/(hdr.Fs/1000),1);
    
    %Hardcode in expected ms value
    if strcmp(filenicks{i}, 'onesec') || strcmp(filenicks{i}, 'newpres');
        tempdat(:,7) = (0:1000:(numel(tempdat(:,6))-1)*1000)';
    elseif strcmp(filenicks{i}, 'tensec')
       tempdat(:,7) = (0:10000:(numel(tempdat(:,6))-1)*10000)';   
    end
    
    tempdat(:,8) = tempdat(:,6) - tempdat(:,7);
    
    %Save in struct
    trigsamples.(filenicks{i}) = tempdat;
    
    for ii = 1:numel(tempdat(:,5))-1
        steps(i,ii) = tempdat(ii+1,5) - tempdat(ii,5);
    end
    %clear tempdat
end

%For figuring out ratio of offset errors WIP
numbers = unique(steps(1,:));
count = hist(steps(1,:), numbers);

%% Plots of difference - recorded vs expected
figure; hold on;

subplot(3,2,1);
plot(trigsamples.onesec(:,8))
title({["Actual - Expected trigger time"]})
xlabel("Duration (sec)");
ylabel("Difference (ms)");

subplot(3,2,2);
plot(trigsamples.onesec(:,5), 'Color', [0.93 0.7 0.12])
title({["Actual - Expected trigger time"]})
xlabel("Duration (sec)");
ylabel("Difference (n samples)");

subplot(3,2,3);
plot(trigsamples.newpres(:,8))
title({["Actual - Expected trigger time"], ["new presentation version"]})
xlabel("Duration (sec)");
ylabel("Difference (ms)");

subplot(3,2,4);
plot(trigsamples.newpres(:,5), 'Color', [0.93 0.7 0.12])
title({["Actual - Expected trigger time"], ["new presentation version"]})
xlabel("Duration (sec)");
ylabel("Difference (n samples)");

subplot(3,2,5);
plot(trigsamples.tensec(:,8))
title({["Actual - Expected trigger time"], ["trigger every 10 sec"]})
xlabel("Duration (10 sec)");
ylabel("Difference (ms)");

subplot(3,2,6);
plot(trigsamples.tensec(:,5), 'Color', [0.93 0.7 0.12])
title({["Actual - Expected trigger time"], ["trigger every 10 sec"]})
xlabel("Duration (10 sec)");
ylabel("Difference (n samples)");
