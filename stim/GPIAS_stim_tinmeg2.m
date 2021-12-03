% Create GPIAS trials in broadband noise
% Author: Niklas Edvall, niklas.edvall@ki.se

%tinmeg2 is updated with options for NBN carrier noise and pure tones to "simulate tinnitus" in silent gaps.
%some code for HP/LP-filter designs are different from tinmeg1 that was mainly ran on Matlab R2019-2020 (local PC Biomedicum).

% soundfile is saved to cd/output as 44.1kHz, 16bit (default) wav as filename

filename = 'BP3_GPi240.wav';

% Variables to specify:
fs = 44100;         % Hz, samplerate (samples per second)
dt = 1/fs;          % seconds per sample

lowpf = 18000;      %Lowpass filter cutoff

useNBN = 1;         %1 to use narrow band noise as carrier - default is white noise
bpfiltfreq = 3000;  %Center frequency of band pass filter

gaptone = 'yes';    %'yes' to fill gap with pure tone - default is silent gap
gaptonef = 3000;    %frequency of gap pure tone
gaptonelvlv = 45;   %level of tone in gap

predur = 0.750;     %sec, pre-duration
fallt = 0.002;      %sec, fall-time before
gapdur = 0.050;     %sec, gap duration
riset = 0.002;      %sec, rise-time after silent gap
ISI = 0.240;        %sec, Inter-stimulus interval (end of risetime to pulse)
pulsedur = 0.020;   %sec, instantaneous rise and fall 
postdur = 4;        %sec, post-duration

bkglvl = 60;        %dB, level of carrier noise
pulselvl = 90;      %dB, level of startle pulse

callvl = 95;        %Calibration (i.e maximum level to be presented)
                    %Magnitude = 1, dB = 0
                    %All other levels relative to this

%Bandpass filter for NBN carrier noise
octfilter = octaveFilter(bpfiltfreq, '1/3 octave','SampleRate', fs);

% Signal created from variables above
pulsediff = db2mag((callvl - pulselvl)*-1); 
bkgdiff = db2mag((callvl - bkglvl)*-1); 
% minus one as the calibration level is the reference at outpath 0dB or magnitude 1,
% i.e callvl is unmodified at magnitude 1 and other amplitudes are lowered by level difference (pulsediff or bkgdiff)
% Calibrate to callvl accordingly.

% Create amplitude-window/envelope
pre = ones(1, round(predur*fs)) .* bkgdiff;        %duration and level of prestim window
post = ones(1, round(postdur*fs)) .* bkgdiff;      %duration and level of poststim window
gap = zeros(1, round(gapdur*fs));                  %duration and level of gap window
pulse = ones(1, round(pulsedur*fs)) .* pulsediff;  %duration and level of pulse window
ISI2 = ones(1, round(ISI*fs)) .* bkgdiff;          %duration and level of ISI

% WIP - If gaptone, create pure tone att gaptonef and replace in silent gap
if gaptone == 'yes';
    pt = sin(2*pi*gaptonef*(0:dt:gapdur));
    %gap = pt;
end

rffreq = 1/(fallt * 2);  %Frequency that has period of 2 fallt
t = (0:dt:fallt);        %vector for rise/fall

%General sine is = Amplitude * sin(2*pi*f*t + phase) + Amp shift
%Fall-window of duration fallt with t samples:
fall = (0.5*bkgdiff) * sin(2*pi*rffreq*t + pi/2) + 0.5*bkgdiff; %fall window

rffreq = 1/(riset * 2);  %Frequency that has period of 2 riset
dt = 1/fs;               %seconds per sample
t = (0:dt:riset);        %vector for rise/fall

%General sine is = Amplitude * sin(2*pi*f*t + phase) + Amp shift
%Rise-window of duration riset with t samples:
rise = (0.5*bkgdiff) * sin(2*pi*rffreq*t + pi/2) + 0.5*bkgdiff; %fall window
rise = flip(rise); %rise window

%if NBN - create filtered sections of carrier separately as pulse is always
%white
if useNBN == 1;
    %Create noise of length pre-duartion + fall time
    pren_falln = rand(1, length(pre) + length(fall));
    pren_falln = (pren_falln - 0.5) * 2;
    pren_falln = octfilter(pren_falln');    %Apply NBN filter, octFilt requires signal in column
    pren_falln = pren_falln' .* [pre fall]; %Multiply by envelope to row vector

    risen_ISI2n = rand(1, length(rise) + length(ISI2));
    risen_ISI2n = (risen_ISI2n - 0.5) * 2;
    risen_ISI2n = octfilter(risen_ISI2n');    %Apply NBN filter, octFilt requires signal in column
    risen_ISI2n = risen_ISI2n' .* [rise ISI2];%Multiply by envelope to row vector

    postn = rand(1, length(post));
    postn = (postn - 0.5) * 2;
    postn = octfilter(postn');    %Apply NBN filter, octFilt requires signal in column
    postn = postn' .* [post];     %Multiply by envelope to row vector

    pulsen = rand(1, length(pulse));
    pulsen = (pulsen - 0.5) * 2;          %1.9 for 95% headroom
    pulsen = lowpass(pulsen, lowpf, fs);    %LP filter of noise
    pulsen = pulsen/max(abs(pulsen(:)));    %Limit to 0 +/- 1 range by dividing signal by max()else LP-filter introduces clipping
    pulsen = pulsen .* [pulse];             %Multiply by envelope
end

%Assemble noise if NBN or create noise and fit in window if white carrier
if useNBN == 1;
    noise = [pren_falln gap risen_ISI2n pulsen postn];
    window = [pre fall gap rise ISI2 pulse post];
else
    window = [pre fall gap rise ISI2 pulse post];
    n = rand(1, length(window));  %Create rand vector of same length as window
    n = (n - 0.5) * 2;            %Shift vector to center on zero

    %Lowpass
    n = lowpass(n, lowpf, fs);     %LP filter of noise
    n = n/max(abs(n(:)));          %Limit to 0 +/- 1 range by dividing signal by max()
                                   %else LP-filter introduces clipping

    %Headroom                                
    n = n .* 0.95;                 %Uncomment for 5% headroom
    noise = n .* window;
end

%Output and graph
%Prints total duration of final signal 'noise'
totdur = length(noise)/fs;
totdur = ['Duration of noise is length ', num2str(round(totdur,4)), ' sec'];
disp(totdur);

%Plot amplitude spectrum of final trial + envelope, mind xlim to inspect
%gap/pulse in detail
figure('Position', [100 100 1000 400]); hold on;
plot(window);
plot(noise);
xlim([0.5*fs 1.5*fs]);
set(gca, 'XTick', [0:fs/10:length(noise)]);
set(gca, 'XTickLabel', [0:1/10:5]);
set(gca, 'XGrid', 'on');

saveas(gcf, ['output/' filename '_timespec.svg']);
close;

if useNBN == 1;
    figure('Position', [100 100 1000 400]);
    pspectrum(postn, fs)
    set(gca, 'XScale', 'log');
    xlim([0.1 20]);
    set(gca, 'XTickLabel', [100 1000 10000]);
    xlabel('Frequency (Hz)');
    title('Power spectrum (carrier, post pulse)');

    %saveas(gcf, ['output/' filename '_freqspec.svg']);
    close;
else
    figure('Position', [100 100 1000 400]);
    pspectrum(noise(length(noise)-length(post):end), fs)
    set(gca, 'XScale', 'log');
    xlim([0.1 20]);
    set(gca, 'XTickLabel', [100 1000 10000]);
    xlabel('Frequency (Hz)');
    title('Power spectrum (carrier, post pulse)');

    %saveas(gcf, ['output/' filename '_freqspec.svg']);
    close;
end
    

%audiowrite(['output/' filename], noise, fs);

%% 15 sec Calibration noise - pulsediff and bkgdiff needed! NB: Levels are set as callvl, pulselvl and bkglvl up top!!

caldur = 15;  %sec, Duration of calibration noise

calwinlvl = ones(1, caldur*fs);                 %15 sec window for calibration at calibration level (i.e magnitude = 1)
calwinpulse = ones(1, caldur*fs) .* pulsediff;  %15 sec window for calibration at pulse level
calwinbkg = ones(1, caldur*fs) .* bkgdiff;      %15 sec window for calibration at carrier noise level

caln = rand(1, length(calwinlvl));      %Create rand vector of same length as caldur
caln = (caln - 0.5) * 2;                %Shift vector to center on zero

%Lowpass filter
caln = lowpass(caln, lowpf, fs);      %LP filter of calibration noise
caln = caln/max(abs(caln(:)));        %Limit to 0 +/- 1 range by dividing signal by max()
                                      %else LP-filter introduces clipping
caln = caln .* 0.95;                  %5% headroom for calibration

pspectrum(caln,fs);

calmax = caln .* calwinlvl;           %15 sec of noise at max level
calpulse = caln .* calwinpulse;       %15 sec of noise at pulse level
calbkg = caln .* calwinbkg;           %15 sec of noise at background level

plot(calmax)
axis([0 length(calpulse) min(calpulse)-0.05 max(calpulse)+0.05]);
hold on
plot(calpulse)
plot(calbkg)
hold off

% audiowrite([outpath 'LP10_Calnoise_max' ' (' num2str(callvl) ')' '.wav'], calmax, fs)
% audiowrite([outpath 'LP10_Calnoise_pulse' ' (' num2str(pulselvl) ')' '.wav'], calpulse, fs)
% audiowrite([outpath 'LP10_Calnoise_background' ' (' num2str(bkglvl) ')' '.wav'], calbkg, fs)