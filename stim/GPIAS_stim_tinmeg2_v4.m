
fs = 44100;
dt = 1/fs;

octfilter = octaveFilter(8000, '1/3 octave','SampleRate', fs, 'FilterOrder', 8);

fallt = 0.1;
riset = 0.1;

rffreq = 1/(fallt * 2);  %Frequency that has period of 2 fallt
t = (0+dt:dt:fallt);        %vector for rise/fall

%General sine is = Amplitude * sin(2*pi*f*t + phase) + Amp shift
%Fall-window of duration fallt with t samples:
fall = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %fall window from 1 to 0

rffreq = 1/(riset * 2);  %Frequency that has period of 2 riset
t = (0+dt:dt:riset);        %vector for rise/fall

%General sine is = Amplitude * sin(2*pi*f*t + phase) + Amp shift
%Rise-window of duration riset with t samples:
rise = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %rise window
rise = flip(rise); %rise window

n = (rand(1, 0.8*fs) - 0.5) * 2;

n = octfilter(n');
n = n';

risen = n(1:length(rise)) .* rise;
falln = n(end-length(fall)+1:end) .* fall;

cf = risen+falln;

n_cut = n(4411:end-4410);

figure; subplot(2,1,1);
spectrogram([cf n_cut cf n_cut cf n_cut cf n_cut cf n_cut], 'yaxis', 800, 120, 600, fs, 'MinThreshold', -110);
title('With crossfade between trials');
subplot(2,1,2);
spectrogram([n_cut n_cut n_cut n_cut n_cut n_cut], 'yaxis', 800, 120, 600, fs, 'MinThreshold', -110);
title('No crossfade between trials')
