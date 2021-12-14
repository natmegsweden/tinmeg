
fs = 44100;          % sample rate
dt = 1/fs;           % seconds per sample

ptf = 1000;          % Carrier tone frequency (Hz)

totdur = 1;          % total duration of stimuli (sec)

modfreq1 = 31;        % Frequency for amplitude modulation 1
moddepth1 = 1;        % Modulation depth of amplitude modulation 1

modfreq2 = 41;        % Frequency for amplitude modulation 1
moddepth2 = 1;        % Modulation depth of amplitude modulation 1

modfreq3 = 47;        % Frequency for amplitude modulation 1
moddepth3 = 1;        % Modulation depth of amplitude modulation 1

riset = 0.1;            % rise time at start of stimuli (sec)
fallt = 0.1;            % fall time at end of stimuli (sec)

modamp1 = 1;      % Amplitude of modulation 1
modamp2 = 1;      % Amplitude of modulation 2
modamp3 = 1;      % Amplitude of modulation 3

t = (0:dt:totdur-dt);   %n samples for total duration

carrieramp = modamp1/moddepth1; % Amplitude of carrier tone

pt = carrieramp * sin(2*pi*ptf*t); % Create carrier tone

%Rise/fall envelope
rffreq = 1/(fallt * 2);  %Frequency that has period of 2 fallt
ft = (0:dt:fallt);        %vector for rise/fall

%Fall/rise-window of duration fallt with t samples:
fall = 0.5 * sin(2*pi*rffreq*ft + pi/2) +0.5; %fall window
rise = flip(fall);
RFenv = [rise ones(1, length(t)-(length(fall)+length(rise))) fall];

ptam = carrieramp*(1+moddepth1*sin(2*pi*modfreq1*t)) .* sin(2*pi*ptf*t);
ptam = carrieramp*(1+moddepth2*sin(2*pi*modfreq2*t)) .* ptam;
ptam = carrieramp*(1+moddepth3*sin(2*pi*modfreq3*t)) .* ptam;
ptam = ptam .* RFenv;

figure; hold on;
plot(ptam);
plot(carrieramp*(1+moddepth2*sin(2*pi*modfreq2*t)))
plot(carrieramp*(1+moddepth1*sin(2*pi*modfreq1*t)))
plot(carrieramp*(1+moddepth3*sin(2*pi*modfreq3*t)))
plot(RFenv)
legend('Signal (1000 Hz)', 'AM envelope 1 (31 Hz)', 'AM envelope 2 (41 Hz)', 'AM envelope 3 (47 Hz)', 'Rise/Fall envelope')