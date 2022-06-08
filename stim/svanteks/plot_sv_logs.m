%% Import data from text file.
filename = 'P:\C3_Canlon lab\Auditory Physiology\Projects\Tinnitus\Human\(9) TinMEG1\Calibration documentation\TinMEG NBN\LOG1395_2msLogstep.csv';
delimiter = {''};

% Format for each line of text:
%   column1: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Create output variable
svlog = [dataArray{1:end-1}];
svlog = svlog';

% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Keep trials of interest, first pulse sample (_start) identified in plot(svlog)

ntrials = 5;

tafter = 100; % 1 sample = 2 ms
tbefore = 169;

%Inspect to find trials
%plot(svlog)

%% 8kHz NBN
GOstart = 15480;
GO_offset = 1000;

GPstart = 20622;
GP_offset = 1500;

POstart = 27976;
PO_offset = 1000;

%GO plot
figure('Position', [280,110,1200,850]); subplot(3,1,1); hold on;
for i = 0:4
    temp(i+1,:) = svlog(146+GOstart-tbefore+GO_offset*i:146+GOstart+tafter+GO_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
%xlabel('Time (ms)');
%ylabel('dB SPL equivalent - 2 ms logstep')
title('Gap, 8kHz NBN + Bob');

%GP plot
subplot(3,1,2); hold on;
for i = 0:4
    temp(i+1,:) = svlog(GPstart-tbefore+GP_offset*i:GPstart+tafter+GP_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
%xlabel('Time (ms)');
ylabel('dB SPL equivalent - 2 ms logstep')
title('Gap + Pulse, 8kHz NBN + Bob');

%PO plot
subplot(3,1,3); hold on;
for i = 0:4
    temp(i+1,:) = svlog(POstart-tbefore+PO_offset*i:POstart+tafter+PO_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
xlabel('Time (ms)');
%ylabel('dB SPL equivalent - 2 ms logstep')
title('Pulse, 8kHz NBN + Bob');

saveas(gcf, 'Bobtrials_NBN8.svg');
close;

%% 3kHz NBN
GOstart = 37017;
GO_offset = 1000;

GPstart = 42162;
GP_offset = 1500;

POstart = 49515;
PO_offset = 1000;

%GO plot
figure('Position', [280,110,1200,850]); subplot(3,1,1); hold on;
for i = 0:4
    temp(i+1,:) = svlog(146+GOstart-tbefore+GO_offset*i:146+GOstart+tafter+GO_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
%xlabel('Time (ms)');
%ylabel('dB SPL equivalent - 2 ms logstep')
title('Gap, 3kHz NBN + Bob');

%GP plot
subplot(3,1,2); hold on;
for i = 0:4
    temp(i+1,:) = svlog(GPstart-tbefore+GP_offset*i:GPstart+tafter+GP_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
%xlabel('Time (ms)');
ylabel('dB SPL equivalent - 2 ms logstep')
title('Gap + Pulse, 3kHz NBN + Bob');

%PO plot
subplot(3,1,3); hold on;
for i = 0:4
    temp(i+1,:) = svlog(POstart-tbefore+PO_offset*i:POstart+tafter+PO_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
xlabel('Time (ms)');
%ylabel('dB SPL equivalent - 2 ms logstep')
title('Pulse, 3kHz NBN + Bob');

saveas(gcf, 'Bobtrials_NBN3.svg');
close;

%% BBN
GOstart = 58557;
GO_offset = 1000;

GPstart = 63702;
GP_offset = 1500;

POstart = 71055;
PO_offset = 1000;

%GO plot
figure('Position', [280,110,1200,850]); subplot(3,1,1); hold on;
for i = 0:4
    temp(i+1,:) = svlog(146+GOstart-tbefore+GO_offset*i:146+GOstart+tafter+GO_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
%xlabel('Time (ms)');
%ylabel('dB SPL equivalent - 2 ms logstep')
title('Gap, BBN + Bob');

%GP plot
subplot(3,1,2); hold on;
for i = 0:4
    temp(i+1,:) = svlog(GPstart-tbefore+GP_offset*i:GPstart+tafter+GP_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
%xlabel('Time (ms)');
ylabel('dB SPL equivalent - 2 ms logstep')
title('Gap + Pulse, BBN + Bob');

%PO plot
subplot(3,1,3); hold on;
for i = 0:4
    temp(i+1,:) = svlog(POstart-tbefore+PO_offset*i:POstart+tafter+PO_offset*i);
    plot(temp(i+1,:), 'Color', [0 0 0 0.25]);
end
plot(mean(temp, 1), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.2);
xlim([0 tafter+tbefore+1]);
ylim([45 90]);
xticks([0:10:tafter+tbefore+1]);
xticklabels([-340:20:tafter*2]);
xlabel('Time (ms)');
%ylabel('dB SPL equivalent - 2 ms logstep')
title('Pulse, BBN + Bob');

saveas(gcf, 'Bobtrials_BBN.svg');
close;
