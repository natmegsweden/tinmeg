% Create GPIAS trials in broadband noise
% Author: Niklas Edvall, niklas.edvall@ki.se

%tinmeg2 is updated with options for NBN carrier noise and pure tones to "simulate tinnitus" in silent gaps.
%some code for HP/LP-filter designs are different from tinmeg1 that was mainly ran on Matlab R2019-2020 (local PC Biomedicum).

% soundfile is saved to cd/output as 44.1kHz, 16bit (default) wav as filename

%% To do

% Print time constants to table
% Automate filename
% Calibration for NBNs

%%
filename = 'BP3_GPi240.wav';

% Variables to specify:
fs = 44100;         % Hz, samplerate (samples per second)
dt = 1/fs;          % seconds per sample

lowpf = 18000;      %Lowpass filter cutoff

useNBN = 0;         %1 to use narrow band noise as carrier - default is white noise
bpfiltfreq = 750;  %Center frequency of band pass filter

gaptone = 1;        %1 to fill gap with pure tone - default is silent gap
gaptonef = 2800;    %frequency of gap pure tone
gaptonelvl = 50;    %level of tone in gap

crossrisefall = 1;  %1 to overlap rise/falltime for carrier and gaptone 50%. Requires fallt and riset to be same

predur = 0.1;     %sec, pre-duration
fallt = 0.002;      %sec, fall-time before
gapdur = 0.055;     %sec, gap duration (or tone if gaptone = 1)
riset = 0.002;      %sec, rise-time after silent gap
ISI = 0.060;        %sec, Inter-stimulus interval (end of risetime to pulse)
pulsedur = 0.020;   %sec, instantaneous rise and fall 
postdur = 0.12;        %sec, post-duration

bkglvl = 65;        %dB, level of carrier noise
pulselvl = 90;      %dB, level of startle pulse

callvl = 95;        %Calibration (i.e maximum level to be presented)
                    %Magnitude = 1, dB = 0
                    %All other levels relative to this

%Bandpass filter for NBN carrier noise
octfilter = octaveFilter(bpfiltfreq, '1/3 octave','SampleRate', fs);

% Signal created from variables above
pulsediff = db2mag((callvl - pulselvl)*-1); 
bkgdiff = db2mag((callvl - bkglvl)*-1); 
gaptonediff = db2mag((callvl - gaptonelvl)*-1);
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
if gaptone == 1;
    
    %pure tone with same duration as specified for gap
    pt = sin(2*pi*gaptonef*(0+dt:dt:gapdur));
    
    %Rise and fall envelopes for puretone (same duration/variable as rise/fall-time
    %for carrier noise
    rffreq = 1/(riset * 2);  %Frequency that has period of 2 riset
    t = (0:dt:riset);        %vector for rise/fall
    ptrise = flip((0.5*gaptonediff) * sin(2*pi*rffreq*t + pi/2) + 0.5*gaptonediff); %PT rise window
    
    rffreq = 1/(fallt * 2);  %Frequency that has period of 2 riset
    t = (0:dt:fallt);        %vector for rise/fall
    ptfall = (0.5*gaptonediff) * sin(2*pi*rffreq*t + pi/2) + 0.5*gaptonediff; %PT fall window
    
    %pure tone envelope
    ptenv = ones(1, round(gapdur*fs)) .* gaptonediff;  %duration and level of gaptone window
    ptenv = ptenv(1:end-(length(ptrise) + length(ptfall)));
    ptenv = [ptrise ptenv ptfall];
        
    %Pad ptenv if dimensions mismatch
    if length(ptenv) > length(pt)
        ptenv = [ptenv(1:length(ptenv)/2) ptenv(length(ptenv)/2+1+(length(ptenv)-length(pt)):end)]
    end
    
    ptn = pt .* ptenv;      %Multiply by envelope
    
end

rffreq = 1/(fallt * 2);  %Frequency that has period of 2 fallt
t = (0:dt:fallt);        %vector for rise/fall

%General sine is = Amplitude * sin(2*pi*f*t + phase) + Amp shift
%Fall-window of duration fallt with t samples:
fall = (0.5*bkgdiff) * sin(2*pi*rffreq*t + pi/2) + 0.5*bkgdiff; %fall window

rffreq = 1/(riset * 2);  %Frequency that has period of 2 riset
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
    pulsen = (pulsen - 0.5) * 2;            
    pulsen = lowpass(pulsen, lowpf, fs);      %LP filter of noise
    pulsen = pulsen/max(abs(pulsen(:)));      %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    pulsen = pulsen .* 0.95;                  % 5% headroom
    pulsen = pulsen .* [pulse];               % Multiply by envelope
    
else %default: use white noise
    pren_falln = rand(1, length(pre) + length(fall));
    pren_falln = (pren_falln - 0.5) * 2;
    pren_falln = lowpass(pren_falln, lowpf, fs);     %LP filter of noise
    pren_falln = pren_falln/max(abs(pren_falln(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    pren_falln = pren_falln .* 0.95;                 % 5% headroom
    pren_falln = pren_falln .* [pre fall];           % Multiply by envelope

    risen_ISI2n = rand(1, length(rise) + length(ISI2));
    risen_ISI2n = (risen_ISI2n - 0.5) * 2;
    risen_ISI2n = lowpass(risen_ISI2n, lowpf, fs);      %LP filter of noise
    risen_ISI2n = risen_ISI2n/max(abs(risen_ISI2n(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    risen_ISI2n = risen_ISI2n .* 0.95;                  % 5% headroom
    risen_ISI2n = risen_ISI2n .* [rise ISI2];           % Multiply by envelope

    postn = rand(1, length(post));
    postn = (postn - 0.5) * 2;
    postn = lowpass(postn, lowpf, fs);      %LP filter of noise
    postn = postn/max(abs(postn(:)));       %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    postn = postn .* 0.95;                  % 5% headroom
    postn = postn .* [post];                % Multiply by envelope

    pulsen = rand(1, length(pulse));
    pulsen = (pulsen - 0.5) * 2;            
    pulsen = lowpass(pulsen, lowpf, fs);      %LP filter of noise
    pulsen = pulsen/max(abs(pulsen(:)));       %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    pulsen = pulsen .* 0.95;                  % 5% headroom
    pulsen = pulsen .* [pulse];                % Multiply by envelope

end

% figure; hold on;
% plot([pre fall ptenv])
% plot([pren_falln ptn])

%Assemble noise
noise = [pren_falln gap risen_ISI2n pulsen postn];
window = [pre fall gap rise ISI2 pulse post];

if gaptone == 1;
    noise = [pren_falln ptn risen_ISI2n pulsen postn];
    window = [pre fall ptenv rise ISI2 pulse post];
    
    %if specified, assemble noise with crossfade of rise/fall (require symmetric rise and fall times)
    if fallt ~= riset
        warning(['Fall and Rise times are different (' num2str(fallt) ' vs ' num2str(riset) ') - not compatible with crossfade of gap']);
    elseif fallt == riset && crossrisefall == 1;
    ol = round(riset/2*fs); %number of samples to overlap, 50% of riset/fallt
    cf1 = pren_falln(end-ol+1:end) + ptn(1:ol); %First overlap region (gap onset, pt start)
    cf2 = ptn(end-ol+1:end) + risen_ISI2n(1:ol); %Second overlap region  (gap offset, pt end)
    
    pren_falln = pren_falln(1:end-ol); %cut the bit getting replaced with cf1
    ptn = ptn(ol+1:end-ol); % +2 ?
    risen_ISI2n = risen_ISI2n(ol+1:end); %cut the bit getting replaced with cf2
    
    noise = [pren_falln cf1 ptn cf2 risen_ISI2n pulsen postn];
    window = [pre fall ptenv rise ISI2 pulse post];
    
    end
    
end


%Output and graph
%Prints total duration of final signal 'noise'
totdur = length(noise)/fs;
totdur = ['Duration of noise is length ', num2str(round(totdur,4)), ' sec'];
disp(totdur);

%Plot amplitude spectrum of final trial + envelope, mind xlim to inspect
%gap/pulse in detail
figure('Position', [100 100 600 400]); hold on;
%plot(window);
plot(noise);
xline(length(pren_falln) + length(cf1)/2);
xline(length(pren_falln) + length(cf1) + length(ptn) + length(cf2)/2)

%xlim([0.7*fs 0.9*fs]);
set(gca, 'XTick', [0:fs/10:length(noise)]);
set(gca, 'XTickLabel', [0:1/10:5]);
set(gca, 'XGrid', 'on');

%saveas(gcf, ['output/' filename '_timespec.svg']);
%close;


figure('Position', [100 100 600 400]);
pspectrum(postn, fs) %noise(1:46042)
set(gca, 'XScale', 'log');
xlim([0.1 20]);
set(gca, 'XTickLabel', [100 1000 10000]);
xlabel('Frequency (Hz)');
title('Power spectrum (carrier, post pulse)');

%saveas(gcf, ['output/' filename '_freqspec.svg']);
%close;

%All variables of final noise in different colors
figure; hold on;
plot(pren_falln);
plot(length(pren_falln)+1:length([pren_falln cf1]), cf1);
plot(length([pren_falln cf1])+1:length([pren_falln cf1 ptn]), ptn);
plot(length([pren_falln cf1 ptn])+1:length([pren_falln cf1 ptn cf2]), cf2);
plot(length([pren_falln cf1 ptn cf2])+1:length([pren_falln cf1 ptn cf2 risen_ISI2n]), risen_ISI2n);
plot(length([pren_falln cf1 ptn cf2 risen_ISI2n])+1:length([pren_falln cf1 ptn cf2 risen_ISI2n pulsen]), pulsen);
plot(length([pren_falln cf1 ptn cf2 risen_ISI2n pulsen])+1:length([pren_falln cf1 ptn cf2 risen_ISI2n pulsen postn]), postn);
xline(length(pren_falln) + length(cf1)/2);


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