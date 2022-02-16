
fs = 44100;

n = (rand(1, 5*fs) - 0.5) * 2;

gap = zeros(1, 0.020*fs);

offset = 0;
for i = 1:5;
    
    a = 1+offset;
    b = fs+offset;
    
    r = round((b-a).*rand(1,1) + a);
    n(r:r+length(gap)-1) = gap;
    
    offset = offset + 44100;
    
end

plot(n);
xticks([1:5*fs/5:5*fs]);
xticklabels([0:5]);
ylim([-1.2 1.2]);
xlim([0 5*fs]);
