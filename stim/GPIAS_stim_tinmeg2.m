% Create GPIAS trials in broadband noise
% Author: Niklas Edvall, niklas.edvall@ki.se

%tinmeg2 is updated with options for NBN carrier noise and pure tones to "simulate tinnitus" in silent gaps.

%some code for HP/LP-filter designs are different from tinmeg1 that was mainly ran on Matlab R2019-2020 (local PC Biomedicum).
%tinmeg2 is made to run on Matlab v9.0 (R2016a) on the NatMEG compute server.

% soundfile is saved to cd/output as 44.1kHz, 16bit (default) wav as filename


%To do:
%recreate lowpass from tinmeg (minimum order stop band 60dB, fpass = 18000, fstop: unknown (check stim)

addpath('stim');

outpath = ['../../stim_output/'];

filename = 'B_C70_P95.wav';

%Load filters


% Variables to specify:
fs = 44100;         % Hz, samplerate (samples per second)
lowpf = 12000;      %Lowpass filter cutoff
dt = 1/fs;          % seconds per sample

gaptone = 'yes';    %Fill gap with pure tone? yes or no
gaptonef = 3000;    %frequency of gap pure tone
gaptonelvlv = 45;   %level of tone in gap

predur = 0.750;     %sec, pre-duration
fallt = 0.052;      %sec, fall-time before
gapdur = 0.050;     %sec, gap duration
riset = 0.002;      %sec, rise-time after silent gap
ISI = 0.120;        %sec, Inter-stimulus interval (end of risetime to pulse)
pulsedur = 0.020;   %sec, instantaneous rise and fall 
postdur = 4;        %sec, post-duration

bkglvl = 60;        %dB, level of carrier noise
pulselvl = 80;      %dB, level of startle pulse

callvl = 95;        %Calibration (i.e maximum level to be presented)
                    %Magnitude = 1, dB = 0
                    %All other levels relative to this

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

%If gaptone, create pure tone att gaptonef and replace in silent gap
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

% n for noise
pren_falln = rand(1, length(pre) + length(fall));
pren_falln = (pren_falln - 0.5) * 2; 

% fileter here
% pren_falln = lowpass(pren_falln, lowpf, fs);

%plot spectrum (test)
[pxx,f] = pspectrum(xTable);

plot(f,pow2db(pxx))
xlim([1000 20000]);
ylim([-80 0]);
grid on

set(gca, 'XScale', 'log');

xlabel('Frequency (Hz)')
set(gca, 'XTick', [1000 10000]);
set(gca, 'XTickLabel', [1000 10000]);

ylabel('Power Spectrum (dB)')
title('Default Frequency Resolution')


pren_falln = pren_falln .* [pre fall];

window = [pre fall gap rise ISI2 pulse post];

n = rand(1, length(window));  %Create rand vector of same length as window
n = (n - 0.5) * 2;            %Shift vector to center on zero

%Lowpass
n = lowpass(n, lowpf, fs);     %LP filter of noise
n = n/max(abs(n(:)));           %Limit to 0 +/- 1 range by dividing signal by max()
                                %else LP-filter introduces clipping

%Headroom                                
n = n .* 0.95;                  %Uncomment for 5% headroom
noise = n .* window;

%pspectrum(n,fs);

% Output and graph

%Prints total duration of final signal 'noise'
totdur = length(noise)/fs;
totdur = ['Duration of noise is length ', num2str(round(totdur,4)), ' sec'];
disp(totdur);

plot(window)
axis([0 length(window) min(window)-0.05 max(window)+0.05]);
hold on

plot(noise)
axis([0 length(noise) min(noise)-0.05 max(noise)+0.05]);
hold off


%audiowrite([outpath 'filnoise2.wav'], filnoise, fs)

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