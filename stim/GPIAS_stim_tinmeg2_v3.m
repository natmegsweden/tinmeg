
%"Tinnitus" conditions
tin = [0 3000 8000];

%Background carrier types
bkg = [0 3000 8000];

%Stimulation types
stim = {'GO', 'PO', 'GP'};

%Specify common parameters
fs = 44100;     % samplerate
dt = 1/fs;      % length of sample
lp = 18000;     % lowpass filter frequency
plvl = 90;      % pulse level
bkglvl = 60;    % bakground carrier level
clvl = 90;      % reference calibration level
gapdur = 0.050; % gap duration
isidur = 0.240; % ISI duration
riset = 0.002;  % rise time
pdur = 0.020;   % pulse duration
fallt = 0.002;  % fall time
tonelvl = 50;   % Tone level (dB)

tonelvldiff = db2mag((clvl - tonelvl)*-1); %Tone level difference relative calibration level
    
%Create 15 sec reference calibration noise. NB!! Same as in function: makenoise
calref = (rand(1, 15*fs) - 0.5) * 2;
calref = lowpass(calref, lp, fs);     %LP filter of noise
calref = calref/max(abs(calref(:)));  %Scale to max or LP may introduce clipping

r = 0; %Count output row for variable table

for i = 1:numel(tin);

    for ii = 1:numel(bkg);
        
        %Create 15 sec calibration-files for each carrier
        calnoise60 =   makenoise(15,   fs,  0,    0,  bkg(ii),  lp,  0,  clvl, bkglvl);
        calnoise70 =   makenoise(15,   fs,  0,    0,  bkg(ii),  lp,  0,  clvl, bkglvl+10);
        calnoise80 =   makenoise(15,   fs,  0,    0,  bkg(ii),  lp,  0,  clvl, bkglvl+20);
        
        audiowrite(['output/audio/bkg_cal60_' num2str(bkg(ii)) '.wav'], calnoise60, fs);
        audiowrite(['output/audio/bkg_cal70_' num2str(bkg(ii)) '.wav'], calnoise70, fs);
        audiowrite(['output/audio/bkg_cal80_' num2str(bkg(ii)) '.wav'], calnoise80, fs);
        clear calnoise60 calnoise70 calnoise80;
        
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
        audiowrite(['output/audio/' padname '.wav'], pad, fs);

        for iii = 1:numel(stim)
            
            r = r + 1; %Increment row counter

            if      strcmp(stim{iii}, 'GO'); rt = riset; ft = fallt;
            elseif  strcmp(stim{iii}, 'PO'); rt = 0; ft = 0; %no rise/fall if PO
            elseif  strcmp(stim{iii}, 'GP'); rt = riset; ft = fallt;
            end

            pre =   makenoise(1,      fs,  0,    ft,  bkg(ii),  lp,  0,  clvl, bkglvl);
            isi =   makenoise(isidur, fs,  rt,   0,   bkg(ii),  lp,  0,  clvl, bkglvl);
            pulse = makenoise(pdur,   fs,  0,    0,   0,        lp,  0,  clvl, plvl);
            post =  makenoise(2.5,    fs,  0,    0,   bkg(ii),  lp,  0,  clvl, bkglvl);
            gap =   zeros(1, fs*gapdur);

            %Assemble noise and specify triggers
            if      strcmp(stim{iii}, 'GO'); n = [pre gap isi post];
                    gapontrig =   length(pre);
                    gapofftrig =  length([pre gap]);
                    pulseontrig = NaN;
            elseif  strcmp(stim{iii}, 'PO'); n = [pre isi pulse post];
                    gapontrig =   NaN;
                    gapofftrig =  NaN;
                    pulseontrig = length([pre isi]);
            elseif  strcmp(stim{iii}, 'GP'); n = [pre gap isi pulse post];
                    gapontrig =   length(pre);
                    gapofftrig =  length([pre gap]);
                    pulseontrig = length([pre gap isi]);
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
            
            subplot(3, 2, [1 2]); hold on
            plot(n);
            if exist('pt', 'var') == 1; plot(pt); end
            title(['Full trial: ' fname], 'Interpreter', 'none');
            xlim([0 length(n)]);
            ylim([-1 1]);
            set(gca, 'XTick', [0:fs/5:length(n)]);
            set(gca, 'XTickLabel', [0:1/5:stimdur]);
            set(gca, 'XGrid', 'on');
            xlabel('Time (sec)');
            if ~isnan(gapontrig); xline(gapontrig); end
            if ~isnan(gapofftrig); xline(gapofftrig); end
            if ~isnan(pulseontrig); xline(pulseontrig); end
            
            if ~strcmp(stim{iii}, 'PO')
                subplot(3, 2, 3); hold on;
                plot(n);
                if exist('pt', 'var') == 1; plot(pt); end
                xlim([length(pre)-(fallt*fs*2) length([pre gap])+fallt*fs*2]);
                set(gca, 'XTick', [0:fs/100:length(n)]);
                set(gca, 'XTickLabel', [0:1/100:stimdur]);
                set(gca, 'XGrid', 'on');
                xlabel('Time (sec)');
                title('Zoom in - gap');
                if ~isnan(gapontrig); xline(gapontrig); end
                if ~isnan(gapofftrig); xline(gapofftrig); end
                if ~isnan(pulseontrig); xline(pulseontrig); end
            end
            
            if ~strcmp(stim{iii}, 'GO')
                subplot(3, 2, 4); hold on;
                plot(n);
                if exist('pt', 'var') == 1; plot(pt); end
                xlim([pulseontrig-(fallt*fs*2) pulseontrig+(length(pulse))+fallt*fs*2]);
                set(gca, 'XTick', [0:fs/100:length(n)]);
                set(gca, 'XTickLabel', [0:1/100:stimdur]);
                set(gca, 'XGrid', 'on');
                xlabel('Time (sec)');
                title('Zoom in - pulse');
                if ~isnan(gapontrig); xline(gapontrig); end
                if ~isnan(gapofftrig); xline(gapofftrig); end
                if ~isnan(pulseontrig); xline(pulseontrig); end
            end
            
            subplot(3, 2, [5 6]);
            pspectrum(n(end-length(post):end), fs);
            set(gca, 'XScale', 'log');
            xlim([0.1 20]);
            set(gca, 'XTickLabel', [100 1000 10000]);
            title('Frequency spectrum - post pulse');
            xlabel('Frequency (Hz)');
            
            %save figures
            saveas(gcf, ['output/figures/' fname '.svg']);
            saveas(gcf, ['output/figures/' fname '.png']);
            close;
            
            %save audiofile
            audiowrite(['output/audio/' fname '.wav'], n, fs);
            
            clear temptin tempbkg fname n pt pre isi pulse post gap
            
        end
    end
end

%Assemble table with variables of interest (i.e. for trigger timings)
varnames = {'Stimulus', 'trial_duration', 'gap_on', 'gap_off', 'pulse_on', 'max_mag'};
stimtab = table(name', stimdur', gapon', gapoff', pulseon', maxmag', 'VariableNames', varnames);

clear i ii iii name stimdur gapon gapoff pulseon maxmag gapofftrig gapontrig r tonelvldiff calref varnames

writetable(stimtab, 'output/stimtable.xlsx');

%% Consider if pre-creating ITI files
%integer number of periods of pt using %fplot(@(x) 1/x/dt, [2800 3200]) and
%nextprime();
