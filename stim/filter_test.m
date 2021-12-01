
fs = 44100;

noise = rand(1, fs*5);
noise = (noise - 0.5) * 1.9;
noise = noise';

thirdoct_CF8 = octaveFilter(8000, '1/3 octave','SampleRate', fs);
thirdoct_CF3 = octaveFilter(3000, '1/3 octave','SampleRate', fs);

filtnoise2 = octFilt(noise);

fs = 44100;
dur = 5;
t = (0:1/Fs:dur)';
x = rand(1, length(t));
x = (x - 0.5) * 1.9;
x = x';

x = octFilt(x);

xTable = timetable(seconds(t),x);

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