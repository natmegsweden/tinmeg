
%To do:

fs = 44100; %Samplerate
dt = 1/fs;  %Seconds per sample

lowpassf = 18000; %Lowpass filter cutoff for calibration noise, pulse and BBN background

ntrials = 5;                %number of presentations per block
tin = [0 3000 8000];        %"Tinnitus" conditions (PT frequencies)
bkg = [0 3000 8000];        %Background carrier types (0: BBN or NBN center frequency)
stim = {'GO', 'PO', 'GP'};  %Stimulation types (Gap Only, Pulse Only, Gap+Pulse)

prepad = 4;         %Duration of background before first trial in each block (sec)
gapdur = 0.050;     %Gap duration (sec)
ISI = 0.240;        %Inter-stimulus-interval, time between gap and pulse in GP trials
pulsedur = 0.020;   %Pulse duration (sec)
totdur = 55;        %Total duration of block, provide as starting point, trimmed later to actual duration (sec)
rf_time = 0.002;    %rise/fall time after/before gap, always symmetric (sec)

minITI = 1.8;       %Minimum inter-trial-interval

pulse_lvl = 90;     %Pulse level (dB)

cal_lvl = 90;       %Reference maximum level, all other are levels relative this.

%Create 15 sec reference calibration noise
calref = (rand(1, 15*fs) - 0.5) * 2;
calref = lowpass(calref, lowpassf, fs); %LP filter of noise
calref = calref/max(abs(calref(:)));    %Scale to max or lowpass may introduce clipping

%Loop for all "tin" conditions
for j = 1:numel(tin)

    %Loop for all background conditions
    for ii = 1:numel(bkg)
        
        bkg_lvl = 60;       %Background noise level (dB)
        tone_lvl = 40;      %Pure tone level in "tin" blocks (dB)
        
        %Name block
        temptin = num2str(tin(j));
        temptin = temptin(1);
        tempbkg = num2str(bkg(ii));
        tempbkg = tempbkg(1);
        fname = ['tin' temptin '_bkg' tempbkg];
        
        %Compensate for equal loudness (ISO 226:2003)
        if bkg(ii) == 8000;
            bkg_lvl = bkg_lvl + 15;
        end
        
        if tin(j) == 8000;
            tone_lvl = tone_lvl + 15;
        end
        
        %Calculate level differences relative calibration level
        bkg_lvldiff = db2mag((cal_lvl - bkg_lvl)*-1);
        pulse_lvldiff = db2mag((cal_lvl - pulse_lvl)*-1);
        tone_lvldiff = db2mag((cal_lvl - tone_lvl)*-1);
        
        %To create calibration file at second level
        cal80diff = db2mag((cal_lvl - bkg_lvl+20)*-1);
    
        %Specify NBN filter parameters
        octfilter = octaveFilter(bkg(ii), '1/3 octave','SampleRate', fs, 'FilterOrder', 8);

        %Construct noise vector to inject trials in
        bkg_noise = (rand(1, totdur*fs) - 0.5) * 2;
        
        %Apply lowpass or NBN filter
        if bkg(ii) > 0;
            bkg_noise = octfilter(bkg_noise'); %Apply NBN filter, octFilt requires signal in column
            bkg_noise = (rms(calref .* bkg_lvldiff)/rms(bkg_noise)) .* bkg_noise; %Scale to match RMS of bakground level reference
            bkg_noise = bkg_noise'; %Pivot back to row vector
        elseif bkg(ii) == 0;
            bkg_noise = lowpass(bkg_noise, lowpassf, fs); %LP filter of noise
            bkg_noise = bkg_noise/max(abs(bkg_noise(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping
            
            bkg_noise80 = (rms(calref .* cal80diff)/rms(bkg_noise)) .* bkg_noise; %For calibration purposes
            
            bkg_noise = (rms(calref .* bkg_lvldiff)/rms(bkg_noise)) .* bkg_noise; %Scale to match RMS of bakground level reference

        end
        
        %Save first 15 sec of bkg_noise from no-tin-blocks for calibration
        if tin(j) == 0 && bkg(ii) == 0
            audiowrite(['output/audio/' fname '_cal60.wav'], bkg_noise(1:15*fs), fs);
            audiowrite(['output/audio/' fname '_cal80.wav'], bkg_noise80(1:15*fs), fs);
            clear bkg_noise80;
        elseif tin(j) == 0 && bkg(ii) > 0
            audiowrite(['output/audio/' fname '_cal60.wav'], bkg_noise(1:15*fs), fs);
        end
        
        %Create rise/fall windows
        rffreq = 1/(rf_time * 2);                     %Frequency that has period of 2 rf_time
        t = (0+dt:dt:rf_time);                        %vector for rise/fall
        fall = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %fall window from 1 to 0
        rise = flip(fall);                            %rise window
        
        %Create and shuffle list of stims to include in block
        stimlist = repmat(stim,1,ntrials);
        Rstimlist = stimlist(randperm(numel(stimlist)));

        offset = prepad*fs; %pad at start of block, offset keep track och trigger times in loop
        r = 1;              %row number for stim/trigger list
        
        %Create an empty vector for trigger times
        PulseOnset = zeros(1, ntrials*numel(stim));
        GapOnset = zeros(1, ntrials*numel(stim));
       
        
        %For each trial
        for i = 1:numel(Rstimlist)
            
            %Create random duration (0-500ms) to add to minimum ITI
            rITI = round(0.5 .* rand(1,1), 2)*fs; %Round to even 10 ms
            
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

        
        tempmatrix = [GOgap', POpulse', GPgap', GPpulse'];
        xlineidx = 0;
        for i = 1:15
            
            trialid = bin2dec(num2str(tempmatrix(i,:)));
            
            if trialid == 8;
                xlineidx = xlineidx + GapOnset(i);
                xline(xlineidx/1000*fs, 'Color', [0 0 0], 'Alpha', 0.3);
            elseif trialid == 4;
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

        %Write stim order and trigger time to table
        varnames = {'Stim', 'STI101_dec', 'GOgap' 'POpulse', 'GPgap', 'GPpulse', 'GO_onset', 'PO_onset'};
        stimtab = table(Rstimlist', STI101_dec', GOgap', POpulse', GPgap', GPpulse', GapOnset', PulseOnset', 'VariableNames', varnames);

        %Write table and soundfile
        writetable(stimtab, ['output/triglists/' fname '.txt'], 'Delimiter', '\t', 'WriteVariableNames', 0);
        audiowrite(['output/audio/' fname '.wav'], stimnoise, fs);

    end
end
