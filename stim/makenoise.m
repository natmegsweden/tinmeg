
function noise = makenoise(duration, samplerate, riset, fallt, nbn_cf, lowpassf, headroom, cal_lvl, output_lvl);

    dt = 1/samplerate;
    lvl_diff = db2mag((cal_lvl - output_lvl)*-1);
    
    if nbn_cf > 0;
        octfilter = octaveFilter(nbn_cf, '1/3 octave','SampleRate', samplerate, 'FilterOrder', 8);
    end
    
    noise = rand(1, duration*samplerate);
    noise = (noise - 0.5) * 2 * (1-headroom);
    
    if nbn_cf > 0;
        noise = octfilter(noise'); %Apply NBN filter, octFilt requires signal in column
        noise = (rms(cal_lvl .* lvl_diff)/rms(noise)) .* noise; %Scale to match RMS of bakground level reference
        noise = noise'; %Pivot back to row vector
    elseif nbn_cf == 0;
        noise = lowpass(noise, lowpassf, samplerate); %LP filter of noise
        noise = noise/max(abs(noise(:))); %Limit to 0 +/- 1 range by dividing signal by max(), else LP-filter introduce clipping                            
        noise = (rms(cal_lvl .* lvl_diff)/rms(noise)) .* noise; %Scale to match RMS of bakground level reference
    end
    
    if riset > 0;
        r_freq = 1/(riset * 2);  %Frequency that has period of 2 riset
        t = (0:dt:riset);        %vector for rise/fall
        
        %Rise-window of duration riset with t samples:
        rise = flip(0.5 * sin(2*pi*r_freq*t + pi/2) + 0.5); %rise window
        
        noise(1:length(rise)) = noise(1:length(rise)) .* rise;
    end
    
    
    if fallt > 0;
        f_freq = 1/(fallt * 2);  %Frequency that has period of 2 fallt
        t = (0:dt:fallt);        %vector for rise/fall
        
        %Fall-window of duration fallt with t samples:
        fall = 0.5 * sin(2*pi*f_freq*t + pi/2) + 0.5; %fall window from 1 to 0
        
        noise(end-length(fall)+1:end) = noise(end-length(fall)+1:end) .* fall;        
    end
  
end

