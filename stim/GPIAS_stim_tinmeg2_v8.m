
% Version 8: Modified from calibration version 5 (i.e. "V5_cal") after
% measuring and calibrating to MSR speaker responses. In V7: set mode to
% create stimulation files or files for calibration purposes (i.e. ordered
% trials with no jitter and no "simulated tinnitus" tone.

% Now include 50ms pre-pulse +10dB over background in broadband background
% only

% Set to 'cal' for calibration or 'exp' for experiment files
mode = 'exp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100; %Samplerate
dt = 1/fs;  %Seconds per sample

lowpassf = 18000; %Lowpass filter cutoff for calibration noise, pulse and BBN background

ntrials = 5;                %number of presentations per block

%"Tinnitus" conditions (PT frequencies)
if mode == 'cal'
    tin = [0];        
elseif mode == 'exp'
    tin = [0 3000 8000];
end

bkg = [0 3000 8000];        %Background carrier types (0: BBN or NBN center frequency)
stim = {'GO', 'PO', 'GP', 'PP'};  %Stimulation types (Gap Only, Pulse Only, Gap+Pulse, Pre-pulse)

prepad = 4;         %Duration of background before first trial in each block (sec)
gapdur = 0.050;     %Gap duration (sec)
ISI = 0.240;        %Inter-stimulus-interval, time between gap and pulse in GP trials
pulsedur = 0.020;   %Pulse duration (sec)
totdur = 90;        %Total duration of block, provide as starting point, trimmed later to actual duration (sec)
rf_time = 0.002;    %rise/fall time after/before gap, always symmetric (sec)

minITI = 2;         %Minimum inter-trial-interval

pulse_lvl = 90-4;   %Pulse level (-2 compensates for MSR speakers)

cal_lvl = 90;       %Reference maximum level, all other are levels relative this.

%Compensate for equal loudness at 8kHz (ISO 226:2003) + 15dB
%AND compensate for speaker frequency response (+3dB at 8kHz, -7 dB at 3kHz)
lvl_comp8 = 15+3;
lvl_comp3 = -7;

%Create 15 sec reference calibration noise
calref = (rand(1, 15*fs) - 0.5) * 2;
calref = lowpass(calref, lowpassf, fs); %LP filter of noise
calref = calref/max(abs(calref(:)));    %Scale to max or lowpass may introduce clipping

%Loop for all "tin" conditions
for j = 1:numel(tin)

    %Loop for all background conditions
    for ii = 1:numel(bkg)
        
        bkg_lvl = 60-7;         %Background noise level (dB), -7 compensates for MSR speakers non-linearity when pulse is amplified to 90dB
        pp_lvl = bkg_lvl+10;    %Pulse level (dB)
        tone_lvl = bkg_lvl-10;  %Pure tone level in "tin" blocks (dB)
        
        %Name block
        temptin = num2str(tin(j));
        temptin = temptin(1);
        tempbkg = num2str(bkg(ii));
        tempbkg = tempbkg(1);
        fname = ['tin' temptin '_bkg' tempbkg];
        
        %Compensate for equal loudness (ISO 226:2003) + 15dB
        %AND compensate for speaker frequency response
        if bkg(ii) == 8000;
            bkg_lvl = bkg_lvl + lvl_comp8;
        elseif bkg(ii) == 3000;
            bkg_lvl = bkg_lvl + lvl_comp3;
        end
        
        if tin(j) == 8000;
            tone_lvl = tone_lvl + lvl_comp8;
        elseif tin(j) == 3000;
            tone_lvl = tone_lvl + lvl_comp3;
        end
        
        %Calculate level differences relative calibration level
        bkg_lvldiff = db2mag((cal_lvl - bkg_lvl)*-1);
        pulse_lvldiff = db2mag((cal_lvl - pulse_lvl)*-1);
        tone_lvldiff = db2mag((cal_lvl - tone_lvl)*-1);
        pp_lvldiff = db2mag((cal_lvl - pp_lvl)*-1);
        
        %To create calibration file at second level
        cal80diff = db2mag((cal_lvl - (bkg_lvl+20))*-1);

        %Construct noise vector to inject trials in
        bkg_noise = (rand(1, totdur*fs) - 0.5) * 2;
        
        %Specify NBN filter parameters
        octfilter = octaveFilter(bkg(ii), '1/3 octave','SampleRate', fs, 'FilterOrder', 8);

        %Apply lowpass or NBN filter
        if bkg(ii) > 0;
            bkg_noise = octfilter(bkg_noise'); %Apply NBN filter, octFilt requires signal in column
            bkg_noise = (rms(calref .* bkg_lvldiff)/rms(bkg_noise)) .* bkg_noise; %Scale to match RMS of bakground level reference
            bkg_noise = bkg_noise'; %Pivot back to row vector
        elseif bkg(ii) == 0;
            bkg_noise = lowpass(bkg_noise, lowpassf, fs); %LP filter of noise
            bkg_noise = bkg_noise/max(abs(bkg_noise(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping
            
            bkg_noise80 = (rms(calref .* cal80diff)/rms(bkg_noise)) .* bkg_noise; %For calibration purposes
            pulse_lvl_cal = (rms(calref .* pulse_lvldiff)/rms(bkg_noise)) .* bkg_noise; %For calibration purposes
            
            bkg_noise = (rms(calref .* bkg_lvldiff)/rms(bkg_noise)) .* bkg_noise; %Scale to match RMS of bakground level reference

        end
        
        %Save first 15 sec of bkg_noise from no-tin-blocks for calibration
        if tin(j) == 0 && bkg(ii) == 0
            audiowrite(['output/audio/' fname '_cal60.wav'], bkg_noise(1:15*fs), fs);
            audiowrite(['output/audio/' fname '_cal80.wav'], bkg_noise80(1:15*fs), fs);
            audiowrite(['output/audio/pulse_lvl_cal.wav'], pulse_lvl_cal(1:15*fs), fs);
            clear bkg_noise80 pulse_lvl_cal;
        elseif tin(j) == 0 && bkg(ii) == 3000
            audiowrite(['output/audio/' fname '_cal60.wav'], bkg_noise(1:15*fs), fs);
        elseif tin(j) == 0 && bkg(ii) == 8000
            audiowrite(['output/audio/' fname '_cal75.wav'], bkg_noise(1:15*fs), fs);
        end
        
        %Create rise/fall windows
        rffreq = 1/(rf_time * 2);                     %Frequency that has period of 2 rf_time
        t = (0+dt:dt:rf_time);                        %vector for rise/fall
        fall = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %fall window from 1 to 0
        rise = flip(fall);                            %rise window
        
        %Create and shuffle list of stims to include in block, pre-pulse
        %only included if bkg and tin == 0.
        if bkg(ii) == 0 && tin(j) == 0;
            stimlist = repmat(stim,1,ntrials);
        else
            stimlist = repmat(stim(1:3),1,ntrials);
        end
        
        %randomize or order stimlist depending on mode
        if mode == 'cal'
            Rstimlist = sort(stimlist);
        elseif mode == 'exp'
            Rstimlist = stimlist(randperm(numel(stimlist)));
        end

        offset = prepad*fs; %pad at start of block, offset keep track of trigger times in loop
        r = 1;              %row number for stim/trigger list
        
        %Create an empty vectors for trigger boolean and wait times
        GOgap = zeros(1, numel(stimlist));
        POpulse = zeros(1, numel(stimlist));
        GPgap = zeros(1, numel(stimlist));
        GPpulse = zeros(1, numel(stimlist));
        PPgap = zeros(1, numel(stimlist));
        PPpulse = zeros(1, numel(stimlist));
        GapOnset = zeros(1, numel(stimlist));
        PulseOnset = zeros(1, numel(stimlist));
        
        %For each trial
        for i = 1:numel(Rstimlist)
            
            %Create random duration (0-500ms) to add to minimum ITI, if
            %mode == cal, no jitter.

            if mode == 'cal'
                rITI = 0; %round(0.5 .* rand(1,1), 2)*fs; %Round to even 10 ms
            elseif mode == 'exp'
                rITI = round(0.5 .* rand(1,1), 2)*fs; %Round to even 10 ms
            end
            
            %Pad offset samples to keep offset at integer ms - avoids cumulative rounding errors in block.
            %https://se.mathworks.com/matlabcentral/answers/440703
            roundpad = mod(-mod(offset,fs),fs);
            
            offsetdiff = offset;

            %if trial is gap only, inject GO trial
            if Rstimlist{i} == 'GO'
                disp(Rstimlist{i});
                %offsetdiff = offset;
                offset = offset + minITI*fs + rITI + roundpad; %update offset
                GO = [fall zeros(1, gapdur*fs) rise];
                bkg_noise(offset:offset+numel(GO)-1) = bkg_noise(offset:offset+numel(GO)-1) .* GO; %inject at offset
                
                %Pad offset samples to keep offset at integer ms - avoids cumulative rounding errors in block.
                %https://se.mathworks.com/matlabcentral/answers/440703
                %disp(['Samples padded: ' num2str(mod(-mod(offset,fs),fs))]);
                %offset = offset+mod(-mod(offset,fs),fs);
                
                if i == 1;
                    offsetdiff = offset;
                elseif i > 1;
                    offsetdiff = offset-offsetdiff;
                end

                %Log trigger time
                GOgap(r) = 1;
                POpulse(r) = 0;
                GPgap(r) = 0;
                GPpulse(r) = 0;
                PPgap(r) = 0;
                PPpulse(r) = 0;
                GapOnset(r) = 1000*(round((offsetdiff/fs),3));
                
            %if trial is pulse only, inject PO trial
            elseif Rstimlist{i} == 'PO'
                disp(Rstimlist{i})
                %offsetdiff = offset;
                offset = offset + minITI*fs + rITI + roundpad; %update offset
                
                %Create a pulse
                PO = (rand(1, pulsedur*fs) - 0.5) * 2;
                PO = lowpass(PO, lowpassf, fs); %LP filter of noise
                PO = PO/max(abs(PO(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
                PO = (rms(calref .* pulse_lvldiff)/rms(PO)) .* PO; %Scale to match RMS of bakground level reference

                bkg_noise(offset:offset+numel(PO)-1) = PO; %Inject at offset

                if i == 1;
                    offsetdiff = offset;
                elseif i > 1;
                    offsetdiff = offset-offsetdiff;
                end
                
                %Log trigger time
                GOgap(r) = 0;
                POpulse(r) = 1;
                GPgap(r) = 0;
                GPpulse(r) = 0;
                PPgap(r) = 0;
                PPpulse(r) = 0;
                PulseOnset(r) = 1000*(round(offsetdiff/fs,3));
                
            %if trial is gap+pulse only, inject GP trial
            elseif Rstimlist{i} == 'GP'
                disp(Rstimlist{i})
                %offsetdiff = offset;
                offset = offset + minITI*fs + rITI + roundpad; %update offset
                GO = [fall zeros(1, gapdur*fs) rise]; %Create gap
                bkg_noise(offset:offset+numel(GO)-1) = bkg_noise(offset:offset+numel(GO)-1) .* GO; %Inject gap at offset
                
                if i == 1;
                    offsetdiff = offset;
                elseif i > 1;
                    offsetdiff = offset-offsetdiff;
                end                
                
                GOgap(r) = 0;
                POpulse(r) = 0;
                GPgap(r) = 1;
                GPpulse(r) = 1;
                PPgap(r) = 0;
                PPpulse(r) = 0;
                GapOnset(r) = 1000*round((offsetdiff/fs),3); %Log trigger time

                offsetdiff = offset;
                offset = offset + numel(GO) + ISI*fs;  %Update offset - risetime is part of ISI
                
                offsetdiff = offset-offsetdiff;
                
                %Create a pulse
                PO = (rand(1, pulsedur*fs) - 0.5) * 2;
                PO = lowpass(PO, lowpassf, fs); %LP filter of noise
                PO = PO/max(abs(PO(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
                PO = (rms(calref .* pulse_lvldiff)/rms(PO)) .* PO; %Scale to match RMS of bakground level reference

                bkg_noise(offset:offset+numel(PO)-1) = PO; %Inject pulse at offset

                PulseOnset(r) = 1000*round(offsetdiff/fs,3); %Log trigger time

             elseif Rstimlist{i} == 'PP'
                disp(Rstimlist{i})
                %offsetdiff = offset;
                offset = offset + minITI*fs + rITI + roundpad; %update offset

                %Create a pre-pulse
                PP = (rand(1, gapdur*fs) - 0.5) * 2;
                PP = lowpass(PO, lowpassf, fs); %LP filter of noise
                PP = PP/max(abs(PP(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
                PP = (rms(calref .* pp_lvldiff)/rms(PP)) .* PP; %Scale to match RMS of bakground level reference

                bkg_noise(offset:offset+numel(PP)-1) = PP; %Inject pp at offset
                
                if i == 1;
                    offsetdiff = offset;
                elseif i > 1;
                    offsetdiff = offset-offsetdiff;
                end
                
                GOgap(r) = 0;
                POpulse(r) = 0;
                GPgap(r) = 0;
                GPpulse(r) = 0;
                PPgap(r) = 1;
                PPpulse(r) = 1;
                GapOnset(r) = 1000*round((offsetdiff/fs),3); %Log trigger time

                offsetdiff = offset;
                offset = offset + numel(PP) + ISI*fs;  %Update offset
                
                offsetdiff = offset-offsetdiff;
                
                %Create a pulse
                PO = (rand(1, pulsedur*fs) - 0.5) * 2;
                PO = lowpass(PO, lowpassf, fs); %LP filter of noise
                PO = PO/max(abs(PO(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
                PO = (rms(calref .* pulse_lvldiff)/rms(PO)) .* PO; %Scale to match RMS of bakground level reference

                bkg_noise(offset:offset+numel(PO)-1) = PO; %Inject pulse at offset

                PulseOnset(r) = 1000*round(offsetdiff/fs,3); %Log trigger time
                
            end
            
            r = r+1;

        end
        
        %Crop full block to size, prepad duration also at end
        stimnoise = bkg_noise(1:offset+prepad*fs);
        
        %If "tinnitus" block, create and add pure tone at tin frequency
        if tin(j) > 0
            pt = sin(2*pi*tin(j)*(0+dt:dt:length(stimnoise)/fs)); %Pure tone of same length as stimnoise
            pt = (rms(calref .* tone_lvldiff)/rms(pt)) .* pt; %Set level of pure tone

            stimnoise = stimnoise + pt;
        end
        
        %Warning if any sample clips
        if max(stimnoise) > 1
            warning(['AMPLITUDE CLIP IN FINAL SOUNDFILE: ' fname])
        end
        
        %Create figure
        figure('units', 'pixels', 'Position', [100 100 800 400]); hold on;
        subplot(3,1,1);
        plot(stimnoise);
        title(['Condition: ' fname], 'interpreter', 'none');
        ylim([-1 1]);
        xlim([0 numel(stimnoise)]);

        set(gca, 'XTick', [0:fs/1:length(stimnoise)]);
        set(gca, 'XTickLabel', [0:1/1:totdur]);

        xlabel('Duration (sec)');
        ylabel('Amplitude');

        tempmatrix = [GOgap', POpulse', GPgap', GPpulse', PPgap', PPpulse'];
        xlineidx = 0;
        for i = 1:numel(Rstimlist)
            
            trialid = bin2dec(num2str(tempmatrix(i,:)));
            
            if trialid == 32;
                xlineidx = xlineidx + GapOnset(i);
                xline(xlineidx/1000*fs, 'Color', [0 0 0], 'Alpha', 0.3);
            elseif trialid == 16;
                xlineidx = xlineidx + PulseOnset(i);
                xline(xlineidx/1000*fs, 'Color', [0 0 0], 'Alpha', 0.3);
            elseif trialid == 12;
                xlineidx = xlineidx + GapOnset(i);
                xline(xlineidx/1000*fs, 'Color', [0 0 0], 'Alpha', 0.3);
                xlineidx = xlineidx + PulseOnset(i);
                xline(xlineidx/1000*fs, 'Color', [0 0 0], 'Alpha', 0.3);
            elseif trialid == 3;
                xlineidx = xlineidx + GapOnset(i);
                xline(xlineidx/1000*fs, 'Color', [0 0 0], 'Alpha', 0.3);
                xlineidx = xlineidx + PulseOnset(i);
                xline(xlineidx/1000*fs, 'Color', [0 0 0], 'Alpha', 0.3);
            end  
        end

        subplot(3,1,2);
        pspectrum(stimnoise(1:prepad*fs), fs)
        ylim([-150 inf]);
        set(gca, 'XScale', 'log');
        xlim([0.1 20]);
        set(gca, 'XTickLabel', [100 1000 10000]);
        title('Frequency spectrum - Pre pad (5 sec)');
        xlabel('Frequency (Hz)');
        
        subplot(3,1,3);
        spectrogram(stimnoise, 'yaxis', 800, 120, 600, fs, 'MinThreshold', -110);
        title('Spectrogram');
        xlabel('Time (sec)');
        
        %save figure
        saveas(gcf, ['output/figures/' fname '.svg']);
        saveas(gcf, ['output/figures/' fname '.png']);
        close;
        
        %Trigger logic        
        if tin(j) == 0
            has_tin = zeros(1, numel(Rstimlist));
        elseif tin(j) > 0
            has_tin = ones(1, numel(Rstimlist));
        end
        
        if tin(j) == 3000
            tin_low = ones(1, numel(Rstimlist));
        elseif tin(j) ~= 3000
            tin_low = zeros(1, numel(Rstimlist));
        end
        
        if bkg(ii) == 0
            has_nbn = zeros(1, numel(Rstimlist));
        elseif bkg(ii) > 0
            has_nbn = ones(1, numel(Rstimlist));
        end
        
        if bkg(ii) == 3000
            nbn_low = ones(1, numel(Rstimlist));
        elseif bkg(ii) ~= 3000
            nbn_low = zeros(1, numel(Rstimlist));
        end
        
        %Decimal trigger padded with zeros in 1s, 2s, 4s and 8s column for 'GOgap' 'POpulse', 'GPgap', 'GPpulse'
        STI101_dec = bin2dec([num2str(has_tin(1)) num2str(tin_low(1)) num2str(has_nbn(1)) num2str(nbn_low(1)), '0', '0', '0', '0']);
        STI101_dec = repmat(STI101_dec,1,numel(Rstimlist));
        
        %Pad last two columns with zero if no PP in conditionm i.e. only if bkg and tin == 0.
        if bkg(ii) ~= 0 && tin(j) ~= 0;
            PPgap = zeros(1, numel(Rstimlist));
            PPgap = zeros(1, numel(Rstimlist));
        end
        
        %varnames = {'Stim', 'STI101_dec', 'GOgap' 'POpulse', 'GPgap', 'GPpulse', 'PPgap', 'PPpulse', 'GO_onset', 'PO_onset'};
        stimtab = table(Rstimlist', STI101_dec', PPgap', PPpulse', GOgap', POpulse', GPgap', GPpulse', GapOnset', PulseOnset');%, 'VariableNames', varnames);

        %Write table and soundfile
        if mode == 'exp'
            %writetable(stimtab, ['output/triglists/' fname '.txt'], 'Delimiter', '\t', 'WriteVariableNames', 0);
            %audiowrite(['output/audio/' fname '.wav'], stimnoise, fs);
        elseif mode == 'cal'
            %writetable(stimtab, ['output/triglists/' fname '_cal.txt'], 'Delimiter', '\t', 'WriteVariableNames', 0);
            %audiowrite(['output/audio/' fname '_cal.wav'], stimnoise, fs);
        end
    end
end
