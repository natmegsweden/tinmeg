% Create GPIAS trials in broadband noise
% Author: Niklas Edvall, niklas.edvall@ki.se

%tinmeg2 is updated with options for NBN carrier noise and pure tones to "simulate tinnitus" in silent gaps.
%some code for HP/LP-filter designs are different from tinmeg1 that was mainly ran on Matlab R2019-2020 (local PC Biomedicum).

% soundfile is saved to cd/output as 44.1kHz, 16bit (default) wav as filename

%% To do

% loop to create files?
% if gaptone == 1; gapdur + riset/2 + fallt/2 <- Needs checking
% Print time constants to table
% Automate filename
% Calibration files

%%
filename = 'BP3_GPi240.wav';

% Variables to specify:
fs = 44100;         % Hz, samplerate (samples per second)
dt = 1/fs;          % seconds per sample

lowpf = 18000;      %Lowpass filter cutoff

useNBN = 1;         %1 to use narrow band noise as carrier - default is white noise
bpfiltfreq = 3000;  %Center frequency of band pass filter

gaptone = 1;        %1 to fill gap with pure tone - default is silent gap
gaptonef = 5000;    %frequency of gap pure tone
gaptonelvl = 50;    %level of tone in gap

crossrisefall = 1;  %1 to overlap rise/falltime for carrier and gaptone 50%.
                    %Requires fallt and riset to be same

predur = 1;       %sec, pre-duration
fallt = 0.002;      %sec, fall-time before
gapdur = 0.050;     %sec, gap duration (or tone if gaptone = 1)
riset = 0.002;      %sec, rise-time after silent gap
ISI = 0.060;        %sec, Inter-stimulus interval (end of risetime to pulse)
pulsedur = 0.020;   %sec, instantaneous rise and fall 
postdur = 1;     %sec, post-duration

bkglvl = 60;        %dB, level of carrier noise
pulselvl = 90;      %dB, level of startle pulse

calrefdur = 5;      %Calibation noise duration to use as reference - WIP
callvl = 90;        %Calibration (i.e maximum level to be presented)
                    %Magnitude = 1, dB = 0
                    %All other levels relative to this

%Only allow crossfade if gaptone is used (or rise/fall times will be shifted)
if gaptone == 0
    crossrisefall = 0;
end

% Bandpass filter for NBN carrier noise (slope in dB/oct = FilterOrder * 6)
octfilter = octaveFilter(bpfiltfreq, '1/3 octave','SampleRate', fs, 'FilterOrder', 8);

% Level differences created from variables above
pulsediff = db2mag((callvl - pulselvl)*-1);
bkgdiff = db2mag((callvl - bkglvl)*-1);
gaptonediff = db2mag((callvl - gaptonelvl)*-1);
% minus one as the calibration level is the reference at output 0dB or magnitude 1,
% i.e callvl is unmodified at magnitude 1 and other amplitudes are lowered by level difference (pulsediff or bkgdiff)
% Calibrate to callvl accordingly.

%Create reference calibration noise
calref = (rand(1, calrefdur*fs) - 0.5) * 2;
calref = lowpass(calref, lowpf, fs);     %LP filter of noise
calref = calref/max(abs(calref(:)));     %Scale to max or LP may introduce clipping

% Create amplitude-window/envelopes
pre = ones(1, round(predur*fs));                   %duration of prestim window
post = ones(1, round(postdur*fs));                 %duration of poststim window
gap = zeros(1, round(gapdur*fs));                  %duration of gap window
pulse = ones(1, round(pulsedur*fs));               %duration of pulse window
ISI2 = ones(1, round(ISI*fs));                     %duration of ISI

%If gaptone, create pure tone att gaptonef and replace in silent gap
if gaptone == 1;
    
    %Compensate duration for 50% overlap if crossfading tone and background noise
    if crossrisefall == 1;
       tonedur = gapdur + (fallt + riset)/2;
    else tonedur = gapdur;
    end
    
    %pure tone with same duration as specified for gap
    pt = sin(2*pi*gaptonef*(0+dt:dt:tonedur));
    
    %Scale to match RMS of gap-tone reference
    pt = (rms(calref .* gaptonediff)/rms(pt)) .* pt;

    %Rise and fall envelopes for puretone (same duration/variable as rise/fall-time
    %for carrier noise
    rffreq = 1/(riset * 2);  %Frequency that has period of 2 riset
    t = (0:dt:riset);        %vector for rise/fall
    ptrise = flip(0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5); %PT rise window
    
    t = (0:dt:fallt);        %vector for rise/fall
    ptfall = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %PT fall window
    
    %pure tone envelope
    ptenv = ones(1, round(tonedur*fs));  %duration and level of gaptone window
    ptenv = ptenv(1:end-(length(ptrise) + length(ptfall))); %duration minus rise and fall
    ptenv = [ptrise ptenv ptfall];

    %Pad ptenv if dimensions mismatch (may occur with short rise/fall due to rounding errors)
    if length(ptenv) > length(pt)
        ptenv = [ptenv(1:length(ptenv)/2) ptenv(length(ptenv)/2+1+(length(ptenv)-length(pt)):end)]
    end
    
    ptn = pt .* ptenv;      %Multiply by envelope
    
end

rffreq = 1/(fallt * 2);  %Frequency that has period of 2 fallt
t = (0:dt:fallt);        %vector for rise/fall

%General sine is = Amplitude * sin(2*pi*f*t + phase) + Amp shift
%Fall-window of duration fallt with t samples:
fall = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %fall window from 1 to 0

rffreq = 1/(riset * 2);  %Frequency that has period of 2 riset
t = (0:dt:riset);        %vector for rise/fall

%General sine is = Amplitude * sin(2*pi*f*t + phase) + Amp shift
%Rise-window of duration riset with t samples:
rise = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %rise window
rise = flip(rise); %rise window

%if NBN - create filtered sections of carrier separately as pulse is always white
if useNBN == 1;

    pren_falln = rand(1, length(pre) + length(fall));       %Create noise of length pre-duartion + fall time
    pren_falln = (pren_falln - 0.5) * 2;                    %Set noise level
    pren_falln = octfilter(pren_falln');                    %Apply NBN filter, octFilt requires signal in column
    pren_falln = (rms(calref .* bkgdiff)/rms(pren_falln)) .* pren_falln; %Scale to match RMS of bakground level reference
    pren_falln = pren_falln';                               %Pivot back to columns
    pren_falln = pren_falln .* [pre fall];                  %Multiply by envelope
    
    risen_ISI2n = rand(1, length(rise) + length(ISI2));     %Create noise of length pre-duartion + fall time
    risen_ISI2n = (risen_ISI2n - 0.5) * 2;                  %Set noise level
    risen_ISI2n = octfilter(risen_ISI2n');                  %Apply NBN filter, octFilt requires signal in column
    risen_ISI2n = (rms(calref .* bkgdiff)/rms(risen_ISI2n)) .* risen_ISI2n; %Scale to match RMS of bakground level reference
    risen_ISI2n = risen_ISI2n';                             %Pivot back to columns
    risen_ISI2n = risen_ISI2n .* [rise ISI2];               %Multiply by envelope

    postn = rand(1, length(post));             %Create noise of length pre-duartion + fall time
    postn = (postn - 0.5) * 2;                 %Set noise level
    postn = octfilter(postn');                 %Apply NBN filter, octFilt requires signal in column
    postn = (rms(calref .* bkgdiff)/rms(postn)) .* postn; %Scale to match RMS of bakground level reference
    postn = postn'; %pivot back to columns
    postn = postn .* post;                     %Multiply by envelope

    pulsen = rand(1, length(pulse));
    pulsen = (pulsen - 0.5) * 2;            
    pulsen = lowpass(pulsen, lowpf, fs);      %LP filter of noise
    pulsen = pulsen/max(abs(pulsen(:)));      %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping
    pulsen = (rms(calref .* pulsediff)/rms(pulsen)) .* pulsen; %Scale to match RMS of pulse level reference;               % Multiply by envelope
    pulsen = pulsen .* pulse;                 %Multiply by envelope

else %default: use white noise
    pren_falln = rand(1, length(pre) + length(fall));
    pren_falln = (pren_falln - 0.5) * 2;
    pren_falln = lowpass(pren_falln, lowpf, fs);     %LP filter of noise
    pren_falln = pren_falln/max(abs(pren_falln(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    pren_falln = (rms(calref .* bkgdiff)/rms(pren_falln)) .* pren_falln; %Scale to match RMS of bakground level reference
    pren_falln = pren_falln .* [pre fall];           %Multiply by envelope

    risen_ISI2n = rand(1, length(rise) + length(ISI2));
    risen_ISI2n = (risen_ISI2n - 0.5) * 2;
    risen_ISI2n = lowpass(risen_ISI2n, lowpf, fs);      %LP filter of noise
    risen_ISI2n = risen_ISI2n/max(abs(risen_ISI2n(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    risen_ISI2n = (rms(calref .* bkgdiff)/rms(risen_ISI2n)) .* risen_ISI2n; %Scale to match RMS of bakground level reference
    risen_ISI2n = risen_ISI2n .* [rise ISI2];           % Multiply by envelope

    postn = rand(1, length(post));
    postn = (postn - 0.5) * 2;
    postn = lowpass(postn, lowpf, fs);      %LP filter of noise
    postn = postn/max(abs(postn(:)));       %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    postn = (rms(calref .* bkgdiff)/rms(postn)) .* postn; %Scale to match RMS of bakground level reference
    postn = postn .* post;                  % Multiply by envelope

    pulsen = rand(1, length(pulse));
    pulsen = (pulsen - 0.5) * 2;            
    pulsen = lowpass(pulsen, lowpf, fs);      %LP filter of noise
    pulsen = pulsen/max(abs(pulsen(:)));      %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    pulsen = (rms(calref .* pulsediff)/rms(pulsen)) .* pulsen; %Scale to match RMS of pulse level reference
    pulsen = pulsen .* pulse;                 % Multiply by envelope

end

%Assemble stimuli
stim = [pren_falln gap risen_ISI2n pulsen postn];

gapontrig = length(pren_falln);
gapofftrig = length([pren_falln gap]);
pulseontrig = length([pren_falln gap risen_ISI2n]);

if gaptone == 1;

    stim = [pren_falln ptn risen_ISI2n pulsen postn];
    
    gapontrig = length(pren_falln);
    gapofftrig = length([pren_falln ptn]);
    pulseontrig = length([pren_falln ptn risen_ISI2n]);

    %if specified, assemble noise with crossfade of rise/fall (require symmetric rise and fall times)
    if fallt ~= riset
        warning(['Fall and Rise times are different (' num2str(fallt) ' vs ' num2str(riset) ') - not compatible with crossfade of gap']);
    elseif fallt == riset && crossrisefall == 1;
    ol = round(riset/2*fs); %number of samples to overlap, 50% of riset/fallt
    cf1 = pren_falln(end-ol+1:end) + ptn(1:ol); %First overlap region (gap onset, pt start)
    cf2 = ptn(end-ol+1:end) + risen_ISI2n(1:ol); %Second overlap region  (gap offset, pt end)
    
    pren_falln = pren_falln(1:end-ol); %cut the bit getting replaced with cf1
    ptn = ptn(ol+1:end-ol); %cut overlapping ends of ptn
    risen_ISI2n = risen_ISI2n(ol+1:end); %cut the bit getting replaced with cf2
    
    stim = [pren_falln cf1 ptn cf2 risen_ISI2n pulsen postn];
    
    gapontrig = length(pren_falln) + (length(cf1)/2);
    gapofftrig = length([pren_falln cf1 ptn]) + (length(cf2)/2);
    pulseontrig = length([pren_falln cf1 ptn cf2 risen_ISI2n]);
    
    end
    
end

%Output and graph
%Prints total duration of final signal 'noise'
totdur = length(stim)/fs;
totdur = ['Duration of stimuli is ', num2str(round(totdur,4)), ' sec'];
disp(totdur);

%% Plot amplitude spectrum of final trial + envelope, mind xlim to inspect

figure('units','normalized','outerposition',[0 0 1 1]); subplot(2,1,1); hold on;
title('stimname');
plot(stim);
xlim([0 length(stim)]);
ylim([-1 1]);
xline(gapontrig);
xline(gapofftrig);
xline(pulseontrig);

set(gca, 'XTick', [0:fs/10:length(stim)]);
set(gca, 'XTickLabel', [0:1/10:5]);
set(gca, 'XGrid', 'on');
xlabel('Time (sec)');

subplot(2,1,2); plot(stim);
xlim([length(pren_falln)-fs/20 length(stim)-length(postn)+fs/20]);
ylim([-1 1]);
xline(gapontrig);
xline(gapofftrig);
xline(pulseontrig);

set(gca, 'XTick', [0:fs/10:length(stim)]);
set(gca, 'XTickLabel', [0:1/10:5]);
set(gca, 'XGrid', 'on');
xlabel('Time (sec)');

%saveas(gcf, ['output/' filename '_timespec.svg']);
%close;

%% Pspectrum fig
figure('units','normalized','outerposition',[0 0 1 1]);
pspectrum(postn, fs);
set(gca, 'XScale', 'log');
xlim([0.1 20]);
set(gca, 'XTickLabel', [100 1000 10000]);
xlabel('Frequency (Hz)');
title(['Power spectrum (carrier, post pulse)' ' --> STIMNAME']);

%saveas(gcf, ['output/' filename '_freqspec.svg']);
%close;

%All variables of final noise in different colors
% figure; hold on;
% plot(pren_falln);
% plot(length(pren_falln)+1:length([pren_falln cf1]), cf1);
% plot(length([pren_falln cf1])+1:length([pren_falln cf1 ptn]), ptn);
% plot(length([pren_falln cf1 ptn])+1:length([pren_falln cf1 ptn cf2]), cf2);
% plot(length([pren_falln cf1 ptn cf2])+1:length([pren_falln cf1 ptn cf2 risen_ISI2n]), risen_ISI2n);
% plot(length([pren_falln cf1 ptn cf2 risen_ISI2n])+1:length([pren_falln cf1 ptn cf2 risen_ISI2n pulsen]), pulsen);
% plot(length([pren_falln cf1 ptn cf2 risen_ISI2n pulsen])+1:length([pren_falln cf1 ptn cf2 risen_ISI2n pulsen postn]), postn);
% xline(length(pren_falln) + length(cf1)/2);

%audiowrite(['output/' filename], noise, fs);

%% 15 sec Calibration noise - pulsediff and bkgdiff needed! NB: Levels are set as callvl, pulselvl and bkglvl up top!!

% rms(pren_falln .* (rms(calbkg)/rms(pren_falln)))

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

%plot(calmax)
axis([0 length(calpulse) min(calpulse)-0.05 max(calpulse)+0.05]);
figure; hold on;
%plot(calpulse)
plot(calbkg)
plot(test)
hold off

% audiowrite([outpath 'LP10_Calnoise_max' ' (' num2str(callvl) ')' '.wav'], calmax, fs)
% audiowrite([outpath 'LP10_Calnoise_pulse' ' (' num2str(pulselvl) ')' '.wav'], calpulse, fs)
% audiowrite([outpath 'LP10_Calnoise_background' ' (' num2str(bkglvl) ')' '.wav'], calbkg, fs)