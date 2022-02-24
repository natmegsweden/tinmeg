
%To do:
%Round onset to ms integer
%ms in table
%tone
%wrap in loop
%save plots for inspection
%clean up/comment

fs = 44100;
dt = 1/fs;

ntrials = 5;                %number of presentations per stim
tin = [0 3000 8000];        %"Tinnitus" conditions
bkg = [0 3000 8000];        %Background carrier types
stim = {'GO', 'PO', 'GP'};  %Stimulation types

prepad = 5;
gapdur = 0.050;
ISI = 0.240;
pulsedur = 0.020;
totdur = 45;

rf_time = 0.002; %rise/fall time

minITI = 1.75; %Minimum ITI

bkg_lvl = 60;
pulse_lvl = 90;

cal_lvl = 90;      % reference maximum level

octfilter = octaveFilter(bkg(ii), '1/3 octave','SampleRate', fs, 'FilterOrder', 8);
lowpassf = 18000;

bkg_lvldiff = db2mag((cal_lvl - bkg_lvl)*-1);
pulse_lvldiff = db2mag((cal_lvl - pulse_lvl)*-1);

%Create 15 sec reference calibration noise. NB!! Same as in function: makenoise
calref = (rand(1, 15*fs) - 0.5) * 2;
calref = lowpass(calref, lowpassf, fs);     %LP filter of noise
calref = calref/max(abs(calref(:)));  %Scale to max or LP may introduce clipping

bkg_noise = (rand(1, totdur*fs) - 0.5) * 2;

if bkg(ii) > 0;
    bkg_noise = octfilter(bkg_noise'); %Apply NBN filter, octFilt requires signal in column
    bkg_noise = (rms(calref .* bkg_lvldiff)/rms(bkg_noise)) .* bkg_noise; %Scale to match RMS of bakground level reference
    bkg_noise = bkg_noise'; %Pivot back to row vector
elseif bkg(ii) == 0;
    bkg_noise = lowpass(bkg_noise, lowpassf, fs); %LP filter of noise
    bkg_noise = bkg_noise/max(abs(bkg_noise(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
    bkg_noise = (rms(calref .* bkg_lvldiff)/rms(bkg_noise)) .* bkg_noise; %Scale to match RMS of bakground level reference
end

rffreq = 1/(rf_time * 2);  %Frequency that has period of 2 fallt
t = (0+dt:dt:rf_time);     %vector for rise/fall

fall = 0.5 * sin(2*pi*rffreq*t + pi/2) + 0.5; %fall window from 1 to 0
rise = flip(fall); %rise window

stimlist = repmat(stim,1,ntrials);
Rstimlist = stimlist(randperm(numel(stimlist)));

offset = 4.75*fs; %pad at start of condition (sec)
r = 1; %row number for stim/trigger list

PulseOnset = zeros(1, ntrials*numel(stim));
GapOnset = zeros(1, ntrials*numel(stim));
GP_G_Onset = zeros(1, ntrials*numel(stim));
GP_P_Onset = zeros(1, ntrials*numel(stim));

for i = 1:numel(Rstimlist)
   
    rITI = floor(round(0.5 .* rand(1,1), 3)*fs); %Rand ITI 0-500 ms, round to integer millisecond
    
    if Rstimlist{i} == 'GO'
        disp(Rstimlist{i});
        offset = offset + minITI*fs + rITI;
        GO = [fall zeros(1, gapdur*fs) rise];
        bkg_noise(offset:offset+numel(GO)-1) = bkg_noise(offset:offset+numel(GO)-1) .* GO;
        
        GapOnset(r) = (offset+rf_time*fs)/fs;
        
    elseif Rstimlist{i} == 'PO'
        disp(Rstimlist{i})
        offset = offset + minITI*fs + rITI;
        
        PO = (rand(1, pulsedur*fs) - 0.5) * 2;
        PO = lowpass(PO, lowpassf, fs); %LP filter of noise
        PO = PO/max(abs(PO(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
        PO = (rms(calref .* pulse_lvldiff)/rms(PO)) .* PO; %Scale to match RMS of bakground level reference
        
        bkg_noise(offset:offset+numel(PO)-1) = PO;
                
        PulseOnset(r) = offset/fs;
        
    elseif Rstimlist{i} == 'GP'
        disp(Rstimlist{i})
        offset = offset + minITI*fs + rITI;
        GO = [fall zeros(1, gapdur*fs) rise];
        bkg_noise(offset:offset+numel(GO)-1) = bkg_noise(offset:offset+numel(GO)-1) .* GO;
        
        GP_G_Onset(r) = (offset+rf_time*fs)/fs;
        
        offset = offset + numel(GO)-(floor(rf_time*fs)) + ISI*fs; %Rise is part of ISI
        
        PO = (rand(1, pulsedur*fs) - 0.5) * 2;
        PO = lowpass(PO, lowpassf, fs); %LP filter of noise
        PO = PO/max(abs(PO(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
        PO = (rms(calref .* pulse_lvldiff)/rms(PO)) .* PO; %Scale to match RMS of bakground level reference
        
        bkg_noise(offset:offset+numel(PO)-1) = PO;
        
        GP_P_Onset(r) = offset/fs;
        
    end
    
    r = r+1;

end



stimnoise = bkg_noise(1:offset+minITI*fs);

figure; hold on;
plot(stimnoise);
for i = 1:numel(GapOnset)
    xline([GapOnset(i)*fs], 'Color', [1 0 0]);
end

for i = 1:numel(PulseOnset)
    xline([PulseOnset(i)*fs], 'Color', [0 1 0]);
end

for i = 1:numel(GP_G_Onset)
    xline([GP_G_Onset(i)*fs], 'Color', [0 0 1]);
end

for i = 1:numel(GP_P_Onset)
    xline([GP_P_Onset(i)*fs], 'Color', [0 0 1]);
end

varnames = {'Stim', 'GO_onset', 'PO_onset', 'GP_Gap_onset', 'GP_Pulse_onset'};
stimtab = table(Rstimlist', GapOnset', PulseOnset', GP_G_Onset', GP_P_Onset', 'VariableNames', varnames);

