% Create sample data (in column vector format) to plot
t = (0:10).';
y1 = t*0.5;
y2 = t;
y3 = t*2;

y = [y1 y2 y3];

load('../mat_data/timelockeds/epochs_eog_all.mat');

% Pack up data
t = (1:165)';
y = epochs_eog_all.PO60{2,5}';
names = regexp(cellstr(sprintf('trial_n_%d ',1:50)),' ','split');
names = names{1,1};

% Plot data
figure('Position',[800 800 1000 800])
hLines = plot(t,y, 'k');
xline([101 101]);
xlim([1 165]);
xticks([1:20:165]);
xticklabels([-500:100:320]);

%save original plot?

% Start brushing mode and wait for user to hit "Enter" when done
brush on
disp('Hit Enter in comand window when done brushing')
pause

% Loop through each graphics object
for k = 1:numel(hLines)
    % Check that the property is valid for that type of object
    % Also check if any points in that object are selected
    if isprop(hLines(k),'BrushData') && any(hLines(k).BrushData)
        % Output the selected data to the base workspace with assigned name
        ptsSelected = logical(hLines(k).BrushData.');
        data = [t(ptsSelected) y(ptsSelected,k)];
        assignin('base',names{k},data)
    end
end

vars = who('trial_n_*');
clear(names{:});

%plot and save trials to be removed in red

close;
clear('colnum');
for i = 1:numel(vars)
    colnum(i) = find(ismember(names, vars{i}));
end

% Plot data
figure('Position',[800 800 1000 1000]); hold on;
plot(t,y, 'k');
plot(t,y(:,[colnum]), 'r');
xline([101 101]);
xlim([1 165]);
xticks([1:20:165]);
xticklabels([-500:100:320]);

y(:,[colnum]) = [];



