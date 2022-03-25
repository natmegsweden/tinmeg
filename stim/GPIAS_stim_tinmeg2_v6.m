
%To do:
%Jitter "mid" och "post" portion of noise for different POGO trials

%"Tinnitus" conditions
tin = [0 3000 8000];

%Background carrier types
bkg = [0 3000 8000];

%Stimulation types
stim = {'GP', 'POG'};

%Specify common parameters
fs = 44100;     % samplerate
dt = 1/fs;      % length of sample
lp = 18000;     % lowpass filter frequency
plvl = 90;      % pulse level

clvl = 90;      % reference maximum level
gapdur = 0.050; % gap duration
isidur = 0.240; % ISI duration
riset = 0.002;  % rise time
pdur = 0.020;   % pulse duration
fallt = 0.002;  % fall time

%Create 15 sec reference calibration noise. NB!! Same as in function: makenoise
calref = (rand(1, 15*fs) - 0.5) * 2;
calref = lowpass(calref, lp, fs);     %LP filter of noise
calref = calref/max(abs(calref(:)));  %Scale to max or LP may introduce clipping

r = 0; %Count output row for variable table

for i = 1:numel(tin);

    for ii = 1:numel(bkg);
        
        bkglvl = 60;    % bakground carrier level
        tonelvl = 40;   % Tone level (dB)
        
        %Compensate for equal loudness (ISO 226:2003)
        if bkg(ii) == 8000;
            bkglvl = bkglvl + 15;
        end
        
        if tin(i) == 8000;
            tonelvl = tonelvl + 15;
        end
        
        tonelvldiff = db2mag((clvl - tonelvl)*-1); %Tone level difference relative calibration level
        
        %Create 15 sec calibration-files (80 only needed for BBN)
        calnoise60 =   makenoise(15,   fs,  0,    0,  bkg(ii),  lp,  0,  clvl, bkglvl);
        %audiowrite(['output/audio/bkg_cal60_' num2str(bkg(ii)) '.wav'], calnoise60, fs);
        
        if bkg(ii) == 0
        calnoise80 =   makenoise(15,   fs,  0,    0,  bkg(ii),  lp,  0,  clvl, bkglvl+20);
        %audiowrite(['output/audio/bkg_cal80_' num2str(bkg(ii)) '.wav'], calnoise80, fs);
        end
        
        clear calnoise60 calnoise80;
        
        %Pad-file (6 sec)
        pad =   makenoise(6,   fs,  0,    0,  bkg(ii),  lp,  0,  clvl, bkglvl);
        
        if tin(i) > 0
            pt = sin(2*pi*tin(i)*(0+dt:dt:length(pad)/fs)); %Pure tone of same length as pad
            pt = (rms(calref .* tonelvldiff)/rms(pt)) .* pt; %Set level of pure tone
            pad = pad + pt;
        end
        
        temptin = num2str(tin(i));
        temptin = temptin(1);
        tempbkg = num2str(bkg(ii));
        tempbkg = tempbkg(1);
        
        padname = ['tin' temptin '_bkg' tempbkg '_pad'];
        
        %Write padfile
        %audiowrite(['output/audio/' padname '.wav'], pad, fs);

        for iii = 1:numel(stim)
            
            r = r + 1; %Increment row counter
            
            pulse = makenoise(pdur,   fs,  0,     0,     0,        lp,  0,  clvl, plvl);
            mid =   makenoise(2.5,    fs,  0,     fallt, bkg(ii),  lp,  0,  clvl, bkglvl);            
            gap =   zeros(1, fs*gapdur);
            isi =   makenoise(isidur, fs,  riset, 0,     bkg(ii),  lp,  0,  clvl, bkglvl);
            post =  makenoise(2.5,    fs,  riset, 0,     bkg(ii),  lp,  0,  clvl, bkglvl);

            %Assemble noise and specify triggers
            if      strcmp(stim{iii}, 'GP'); n = [pulse mid gap isi];
                    pulseontrig = 1;
                    gapontrig =   length([pulse mid]);
            elseif  strcmp(stim{iii}, 'POG'); n = [pulse mid gap post];
                    pulseontrig = 1;
                    gapontrig =   length([pulse mid]);
            end
            
            %Add overlayed tone to noise if "Tinnitus" frequency > 0
            if tin(i) > 0
                pt = sin(2*pi*tin(i)*(0+dt:dt:length(n)/fs)); %Pure tone of same length as noise
                pt = (rms(calref .* tonelvldiff)/rms(pt)) .* pt; %Set level of pure tone
                
                n = n + pt;
            end
            
            temptin = num2str(tin(i));
            temptin = temptin(1);
            tempbkg = num2str(bkg(ii));
            tempbkg = tempbkg(1);
            
            fname = ['tin' temptin '_bkg' tempbkg '_' stim{iii}];
            disp(fname);
            
            %Gather variables for overview in table
            name{1, r} = fname;
            stimdur(r) = length(n)/fs;
            gapon(r) = gapontrig/fs;
            gapoff(r) = gapofftrig/fs;
            pulseon(r) = pulseontrig/fs;
            maxmag(r) = max(n);
                                    
            %% output figure
            figure('units','normalized','outerposition',[0 0 1 1]);
            
            subplot(3, 1, 1); hold on
            plot(n);
            if exist('pt', 'var') == 1; plot(pt); end
            title(['Full trial: ' fname], 'Interpreter', 'none');
            xlim([0 length(n)]);
            ylim([-1 1]);
            set(gca, 'XTick', [0:fs/5:length(n)]);
            set(gca, 'XTickLabel', [0:1/5:stimdur(r)]);
            set(gca, 'XGrid', 'on');
            xlabel('Time (sec)');
            if ~isnan(gapontrig); xline(gapontrig); end
            if ~isnan(pulseontrig); xline(pulseontrig); end
            
            subplot(3, 1, 2); hold on;
            plot(n)
            title(['Gap zoom'], 'Interpreter', 'none');
            xlim([length([pulse mid])-gapdur*fs length([pulse mid gap])+gapdur*fs]);
            ylim([-1 1]);
            set(gca, 'XTick', [0:fs/5:length(n)]);
            set(gca, 'XTickLabel', [0:1/5:stimdur(r)]);
            set(gca, 'XGrid', 'on');
            xlabel('Time (sec)');
            if ~isnan(gapontrig); xline(gapontrig); end
            
            subplot(3, 1, 3);
            pspectrum(n(length(pulse)+1:length(pulse)+length(mid)), fs);
            set(gca, 'XScale', 'log');
            xlim([0.1 20]);
            set(gca, 'XTickLabel', [100 1000 10000]);
            title('Frequency spectrum - post-pulse, pre-gap');
            xlabel('Frequency (Hz)');
            
            %save figures
            %saveas(gcf, ['output/figures/' fname '.svg']);
            %saveas(gcf, ['output/figures/' fname '.png']);
            %close;
            
            %save audiofile
            %audiowrite(['output/audio/' fname '.wav'], n, fs);
            
            %clear temptin tempbkg fname n pt pre isi pulse post gap
            
        end
    end
end

%Assemble table with variables of interest (i.e. for trigger timings)
varnames = {'Stimulus', 'trial_duration', 'gap_on', 'gap_off', 'pulse_on', 'max_mag'};
stimtab = table(name', stimdur', gapon', gapoff', pulseon', maxmag', 'VariableNames', varnames);

%clear i ii iii name stimdur gapon gapoff pulseon maxmag gapofftrig gapontrig r tonelvldiff calref varnames

%writetable(stimtab, 'output/stimtable.xlsx');

%% Consider if pre-creating ITI files
%integer number of periods of pt using %fplot(@(x) 1/x/dt, [2800 3200]) and
%nextprime();
