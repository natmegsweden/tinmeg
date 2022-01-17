% Create GPIAS trials in broadband noise
% Author: Niklas Edvall, niklas.edvall@ki.se

%tinmeg2 is updated with options for NBN carrier noise and pure tones to "simulate tinnitus" in silent gaps.
%some code for HP/LP-filter designs are different from tinmeg1 that was mainly ran on Matlab R2019-2020 (local PC Biomedicum).

% soundfile is saved to cd/output as 44.1kHz, 16bit (default) wav as filename

%% To do

% Same pren length for all stim - simple cut?
% ISI length adapt to triggerdiff?

% Tone level?
% Calibration files

%% Loop to create TinMEG2 stimuli

%Stimulus variables/labels to loop
bkg_type = {'BBN', 'NB3', 'NB8'}; %Carrier noise types
stim_type = {'GO', 'PO', 'GP', 'T3P', 'T8P', 'T3', 'T8'}; %Stimulus type (gap only, pulse only, gap + pulse, tone+pulse)

for i = 1:numel(bkg_type);

    for ii = 1:numel(stim_type);

    filename = [bkg_type{i} '_' stim_type{ii} '.wav'];

    disp(filename)
    
    % Variables to specify:
    fs = 44100;         % Hz, samplerate (samples per second)
    dt = 1/fs;          % seconds per sample

    lowpf = 18000;      %Lowpass filter cutoff

    useNBN = 0;         %1 to use narrow band noise as carrier - default is white noise
    bpfiltfreq = 3000;  %Center frequency of band pass filter
    
    if strcmp(bkg_type{i}, 'NB3');
        useNBN = 1;
        bpfiltfreq = 3000;
    elseif strcmp(bkg_type{i}, 'NB8');
        useNBN = 1;
        bpfiltfreq = 8000;  
    end
    
    gaptone = 0;        %1 to fill gap with pure tone - default is silent gap
    gaptonef = 3000;    %frequency of gap pure tone
    gaptonelvl = 50;    %level of tone in gap

    crossrisefall = 0;  %1 to overlap rise/falltime for carrier and gaptone 50%.
                        %Requires fallt and riset to be same
    
    if strcmp(stim_type{ii}, 'T3') | strcmp(stim_type{ii}, 'T3P')
        gaptone = 1;        %1 to fill gap with pure tone - default is silent gap
        gaptonef = 3000;    %frequency of gap pure tone
        gaptonelvl = 50;    %level of tone in gap

        crossrisefall = 1;  %1 to overlap rise/falltime for carrier and gaptone 50%.
                            %Requires fallt and riset to be same
    elseif strcmp(stim_type{ii}, 'T8') | strcmp(stim_type{ii}, 'T8P')
        gaptone = 1;        %1 to fill gap with pure tone - default is silent gap
        gaptonef = 8000;    %frequency of gap pure tone
        gaptonelvl = 50;    %level of tone in gap

        crossrisefall = 1;  %1 to overlap rise/falltime for carrier and gaptone 50%.
                            %Requires fallt and riset to be same
    end

    predur = 1;         %sec, pre-duration
    fallt = 0.002;      %sec, fall-time before
    gapdur = 0.050;     %sec, gap duration (or tone if gaptone = 1)
    riset = 0.002;      %sec, rise-time after silent gap
    ISI = 0.060;        %sec, Inter-stimulus interval (end of risetime to pulse)
    pulsedur = 0.020;   %sec, instantaneous rise and fall 
    postdur = 1;        %sec, post-duration

    bkglvl = 75;        %dB, level of carrier noise
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

    %If gaptone, create pure tone att gaptonef to replace in silent gap
    if gaptone == 1;
        
        tonedur = gapdur;
        
        %if crossrisefall, compensate for crossover
        if crossrisefall == 1
            tonedur = gapdur + (riset + fallt)/4; 
        end

        %pure tone with tonedur (compensated for crossfade)
        pt = sin(2*pi*gaptonef*(0+dt:dt:tonedur));

        %Scale to match RMS of gap-tone reference
        pt = (rms(calref .* gaptonediff)/rms(pt)) .* pt;
        
        %Rise and fall envelopes for puretone
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
        
        %Split to simplify assembly (i.e of GO trials)
        pren = pren_falln(1:length(pre));
        falln = pren_falln(end-length(fall)+1:end);

        risen_ISI2n = rand(1, length(rise) + length(ISI2));     %Create noise of length pre-duartion + fall time
        risen_ISI2n = (risen_ISI2n - 0.5) * 2;                  %Set noise level
        risen_ISI2n = octfilter(risen_ISI2n');                  %Apply NBN filter, octFilt requires signal in column
        risen_ISI2n = (rms(calref .* bkgdiff)/rms(risen_ISI2n)) .* risen_ISI2n; %Scale to match RMS of bakground level reference
        risen_ISI2n = risen_ISI2n';                             %Pivot back to columns
        risen_ISI2n = risen_ISI2n .* [rise ISI2];               %Multiply by envelope
        
        %Split to simplify assembly (i.e of GO trials)
        risen = risen_ISI2n(1:length(rise));
        ISI2n = risen_ISI2n(end-length(ISI2)+1:end);

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
        
        %Split to simplify assembly (i.e of GO trials)
        pren = pren_falln(1:length(pre));
        falln = pren_falln(end-length(fall)+1:end);

        risen_ISI2n = rand(1, length(rise) + length(ISI2));
        risen_ISI2n = (risen_ISI2n - 0.5) * 2;
        risen_ISI2n = lowpass(risen_ISI2n, lowpf, fs);      %LP filter of noise
        risen_ISI2n = risen_ISI2n/max(abs(risen_ISI2n(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
        risen_ISI2n = (rms(calref .* bkgdiff)/rms(risen_ISI2n)) .* risen_ISI2n; %Scale to match RMS of bakground level reference
        risen_ISI2n = risen_ISI2n .* [rise ISI2];           % Multiply by envelope

        %Split to simplify assembly (i.e of GO trials)
        risen = risen_ISI2n(1:length(rise));
        ISI2n = risen_ISI2n(end-length(ISI2)+1:end);
        
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
    if gaptone == 0;
        stim = [pren falln gap risen ISI2n pulsen postn];

        gapontrig = length([pren falln]);
        gapofftrig = length([pren falln gap]);
        pulseontrig = length([pren_falln gap risen ISI2n]);
    end
    
    if strcmp(stim_type{ii}, 'PO')
        stim = [pren pulsen postn];
        
        gapontrig = NaN;
        gapofftrig = NaN;
        pulseontrig = length([pren]);
    end
    
    if strcmp(stim_type{ii}, 'GO')
        stim = [pren falln gap risen postn];
        
        gapontrig = length([pren falln]);
        gapofftrig = length([pren falln gap]);
        pulseontrig = NaN;
    end

    %if specified, assemble noise with crossfade of rise/fall (require symmetric rise and fall times)
    if fallt ~= riset
        warning(['Fall and Rise times are different (' num2str(fallt) ' vs ' num2str(riset) ') - not compatible with crossfade of gap-tone']);
    elseif fallt == riset && crossrisefall == 1;
        ol = round(riset/2*fs); %number of samples to overlap, 50% of riset/fallt
        
        %figure; hold on; plot(falln); plot(length(falln)+1:length([falln ptn(1:44)]+1), ptn(1:44));
        
        cf1 = falln(end-ol+1:end) + ptn(1:ol); %First overlap region (gap onset, pt start)
        cf2 = ptn(end-ol+1:end) + risen(1:ol); %Second overlap region  (gap offset, pt end)
        
        %figure; hold on; plot(falln(end-ol+1:end)); plot(ptn(1:ol)); plot(cf1);
        
        ptn = ptn(ol+1:end-ol); %cut overlapping ends of ptn
        falln = falln(1:end-ol); %cut overlapping end of pren
        risen = risen(ol+1:end); %cut overlapping beginning of ISI2n

        stim = [pren falln cf1 ptn cf2 risen ISI2n pulsen postn];

        gapontrig = length([pren falln]) + length(cf1)/2;
        gapofftrig = length([pren falln cf1 ptn]) + length(cf2)/2;
        pulseontrig = length([pren falln cf1 ptn cf2 risen ISI2n]);
        
        if strcmp(stim_type{ii}, 'T3') | strcmp(stim_type{ii}, 'T8');
            stim = [pren falln cf1 ptn cf2 risen ISI2n postn];

            gapontrig = length([pren falln]) + length(cf1)/2;
            gapofftrig = length([pren falln cf1 ptn]) + length(cf2)/2;
            pulseontrig = NaN;
        end
        
    end

    %% Plot amplitude spectrum of final trial + envelope, mind xlim to inspect

    figure('units','normalized','outerposition',[0 0 1 1]); subplot(3,1,1); hold on;
    title(filename, 'Interpreter', 'none');
    plot(stim);
    xlim([0 length(stim)]);
    ylim([-1 1]);
    if ~isnan(gapontrig)
        xline(gapontrig); end
    if ~isnan(gapofftrig)
        xline(gapofftrig); end
    if ~isnan(pulseontrig)
        xline(pulseontrig); end

    set(gca, 'XTick', [0:fs/10:length(stim)]);
    set(gca, 'XTickLabel', [0:1/10:5]);
    set(gca, 'XGrid', 'on');
    xlabel('Time (sec)');

    subplot(3,1,2); plot(stim);
    xlim([length([pren falln])-fs/100 length(stim)-length(postn)+fs/100]);
    ylim([-1 1]);
    if ~isnan(gapontrig)
        xline(gapontrig); end
    if ~isnan(gapofftrig)
        xline(gapofftrig); end
    if ~isnan(pulseontrig)
        xline(pulseontrig); end

    set(gca, 'XTick', [0:fs/100:length(stim)]);
    set(gca, 'XTickLabel', [0:1/100:5]);
    set(gca, 'XGrid', 'on');
    xlabel('Time (sec)');

    %Pspectrum fig
    subplot(3,1,3);
    pspectrum(postn, fs);
    set(gca, 'XScale', 'log');
    xlim([0.1 20]);
    set(gca, 'XTickLabel', [100 1000 10000]);
    xlabel('Frequency (Hz)');
    title(['Power spectrum (carrier, post pulse)' filename], 'Interpreter', 'none');

    saveas(gcf, ['output/' filename(1:end-4) '.svg']);
    close;

    %audiowrite(['output/' filename], stim, fs);
    
    %% Output table
    
    %Row number
    if i == 1;
        r = ii;
    elseif i == 2;
        r = ii +(length(stim_type));
    elseif i == 3;
        r = ii +(length(stim_type)*2);
    end;
    
    totdur(r,1) = length(stim)/fs;
    stimname{r,1} = filename;
    gaponlat(r,1) = gapontrig/fs;
    gapofflat(r,1) = gapofftrig/fs;
    pulseonlat(r,1) = pulseontrig/fs;
    gaptrigdiff(r,1) = (gapofftrig - gapontrig)/fs;
    pulsegapdiff(r,1) = (pulseontrig - gapofftrig)/fs;
    
    end
end

%Assmeble table of trigger times
stimtable = table(stimname, totdur, gaponlat, gapofflat, gaptrigdiff, pulseonlat, pulsegapdiff);
writetable(stimtable, 'output/stimtable.xlsx');
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
