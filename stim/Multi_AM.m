
fs = 44100;          % sample rate
dt = 1/fs;           % seconds per sample

ptf = 1000;          % Carrier tone frequency (Hz)

totdur = 1;          % total duration of stimuli (sec)

modfreq1 = 2;        % Frequency for amplitude modulation 1
moddepth1 = 0.2;        % Modulation depth of amplitude modulation 1

modfreq2 = 4;        % Frequency for amplitude modulation 2
moddepth2 = 0.2;        % Modulation depth of amplitude modulation 2

modfreq3 = 8;        % Frequency for amplitude modulation 3
moddepth3 = 0.2;        % Modulation depth of amplitude modulation 3

modfreq4 = 16;        % Frequency for amplitude modulation 4
moddepth4 = 0.2;        % Modulation depth of amplitude modulation 4

modfreq5 = 32;        % Frequency for amplitude modulation 5
moddepth5 = 0.2;        % Modulation depth of amplitude modulation 5

riset = 0.1;            % rise time at start of stimuli (sec)
fallt = 0.1;            % fall time at end of stimuli (sec)

modamp1 = 1;      % Amplitude of modulation 1
modamp2 = 1;      % Amplitude of modulation 2
modamp3 = 1;      % Amplitude of modulation 3
modamp4 = 1;      % Amplitude of modulation 4
modamp5 = 1;      % Amplitude of modulation 5

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
ptam = carrieramp*(1+moddepth4*sin(2*pi*modfreq4*t)) .* ptam;
ptam = carrieramp*(1+moddepth5*sin(2*pi*modfreq5*t)) .* ptam;
ptam = ptam .* RFenv;

%Fix amp-factor for plot (i.e. now fixed at 2000)

figure('Position', [100 100 1200 500]); hold on;
plot(ptam);
plot(2000*sin(2*pi*modfreq1*t))
plot(2000*sin(2*pi*modfreq2*t))
plot(2000*sin(2*pi*modfreq3*t))
plot(2000*sin(2*pi*modfreq4*t))
plot(2000*sin(2*pi*modfreq5*t))
plot(2000*RFenv)
legend('Signal (1000 Hz)', 'AM envelope 1 (2 Hz)', 'AM envelope 2 (4 Hz)', 'AM envelope 3 (8 Hz)', 'AM envelope 4 (16 Hz)', 'AM envelope 5 (32 Hz)', 'Rise/Fall envelope')

%Limit max output to -1/+1
ptam = ptam/max(abs(ptam(:))); 