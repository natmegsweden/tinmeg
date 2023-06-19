%% Sensor level figures och timelockeds

%Load if not already in workspace
if exist('tlk_all_sub','var') == 0
tlk_all_sub = load('../processed_data/timelockeds/aggregated/tlk_all_sub_cmb.mat');
tlk_all_sub = tlk_all_sub.tlk_all_sub;
end

%Load sensor information
load('../input_vars/sensors.mat');

%% Find gradiometers with highest response in PO60_90 TOI

%Prepare structures
L_topgrads = struct();
R_topgrads = struct();

%Define TOI 1 (50 - 150 ms)
toi1 = find(tlk_all_sub.tinmeg1.time == 0.050);
toi2 = find(tlk_all_sub.tinmeg1.time == 0.150);

%Get labels for Right and Left sensors
Rchan_lab = tlk_all_sub.tinmeg1.label(sensors.chanpos(:,1) > 0, :);
Lchan_lab = tlk_all_sub.tinmeg1.label(sensors.chanpos(:,1) < 0, :);

for i = 1:numel(tlk_all_sub.tinmeg1.ID)

    %Separate data in from Right and Left sensors for subject
    Rchan = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(ismember(tlk_all_sub.tinmeg1.label, Rchan_lab),:);
    Lchan = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(ismember(tlk_all_sub.tinmeg1.label, Lchan_lab),:);
    
    %Find sensor on RIGHT side with biggest amplitude response in TOI
    [datR, indR] = sort(mean(Rchan(:,toi1:toi2), 2), 'descend');
    topRname = Rchan_lab{indR(1)};
    topRind = find(ismember(tlk_all_sub.tinmeg1.label, topRname));
    
    %Find sensor on LEFT side with biggest amplitude response in TOI
    [datL, indL] = sort(mean(Lchan(:,toi1:toi2), 2), 'descend');
    topLname = Lchan_lab{indL(1)};
    topLind = find(ismember(tlk_all_sub.tinmeg1.label, topLname));

    %Write top channel index (out of 204)
    L_topgrads.chind{i,:} = topLind;
    R_topgrads.chind{i,:} = topRind;

    %Write top channel label/name
    L_topgrads.chan{i,:} = topLname;
    R_topgrads.chan{i,:} = topRname;

%For subject
end

clear topLname topRname topLind topRind;

%Count up what sensors are most common as top response
[countL, nameL] = groupcounts(L_topgrads.chan);
[countR, nameR] = groupcounts(R_topgrads.chan);

%Find and print info on most common top gradiometer
[C, I] = sort(countL, 'descend');
topL = nameL{I(1)};
disp(['Top sensor LEFT in PO60_90 is: ' topL ' - highest response in ' num2str(countL(I(1))) ' of ' num2str(numel(tlk_all_sub.tinmeg1.ID)) ' participants']);

[C, I] = sort(countR, 'descend');
topR = nameR{I(1)};
disp(['Top sensor RIGHT in PO60_90 is: ' topR ' - highest response in ' num2str(countR(I(1))) ' of ' num2str(numel(tlk_all_sub.tinmeg1.ID)) ' participants']);

clear C I countL countR datL datR L_topgrads R_topgrads Rchan Rchan_lab Lchan Lchan_lab indL indR nameL nameR

%Save topL and topR?

%% PO60_90 tinmeg1

%%%%%%%%%%%%%%%%%%%%%%%
%Average of all grads for each subject
figure; hold on;
for i = 1:numel(tlk_all_sub.tinmeg1.ID)
    plot(mean(tlk_all_sub.tinmeg1.PO60.PO60_90{i}(1:102,:),1), 'Color', [0 0 0 0.45])
    temp_mean(i,:) = mean(tlk_all_sub.tinmeg1.PO60.PO60_90{i}(1:102,:),1);
end
plot(mean(temp_mean), 'Color', [1 0 0 1], 'LineWidth', 2); clear temp_mean

%%%%%%%%%%%%%%%%%%%%%%%
%Top grad L/R
topLind = find(ismember(tlk_all_sub.tinmeg1.label, topL));
topRind = find(ismember(tlk_all_sub.tinmeg1.label, topR));

for i = 1:numel(tlk_all_sub.tinmeg1.ID)
    temp_mean_L(i,:) = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(topLind,:);
    temp_mean_R(i,:) = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(topRind,:);
end

figure; hold on;
plot(mean(temp_mean_L), 'Color', [0 0 1 1], 'LineWidth', 2); clear temp_meanL
plot(mean(temp_mean_R), 'Color', [1 0 0 1], 'LineWidth', 2); clear temp_meanR

%%%%%%%%%%%%%%%%%%%%%%%
%Inhibition in top grad
topLind = find(ismember(tlk_all_sub.tinmeg1.label, topL));
topRind = find(ismember(tlk_all_sub.tinmeg1.label, topR));

for i = 1:numel(tlk_all_sub.tinmeg1.ID)
    temp_PO_L(i,:) = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(topLind,:);
    temp_PO_R(i,:) = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(topRind,:);

    temp_GP_L(i,:) = tlk_all_sub.tinmeg1.GP60.GP60_i240{i}(topLind,:);
    temp_GP_R(i,:) = tlk_all_sub.tinmeg1.GP60.GP60_i240{i}(topRind,:);
end

figure; hold on;
plot(mean(temp_PO_L), 'Color', [0 0 1 1], 'LineWidth', 2); %clear temp_PO_L temp_GP_L
plot(mean(temp_GP_L), '--', 'Color', [0 0 1 1], 'LineWidth', 2); %clear temp_PO_L temp_GP_L

plot(mean(temp_PO_R), 'Color', [1 0 0 1], 'LineWidth', 2); %clear temp_PO_R temp_GP_R
plot(mean(temp_GP_R), '--', 'Color', [1 0 0 1], 'LineWidth', 2); %clear temp_PO_R temp_GP_R

clear temp_GP_R temp_GP_L temp_mean_R temp_mean_L temp_PO_L temp_PO_R

%% tinmeg3

%Inhibition in top grad (!! same sensors as for tinmeg1)
topLind = find(ismember(tlk_all_sub.tinmeg3.label, topL));
topRind = find(ismember(tlk_all_sub.tinmeg3.label, topR));

for i = 1:numel(tlk_all_sub.tinmeg3.ID)
    temp3_PO_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.PO_00{i}(topLind,:);
    temp3_PO_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.PO_00{i}(topRind,:);

    temp3_GP_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.GPP_00{i}(topLind,:);
    temp3_GP_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.GPP_00{i}(topRind,:);
end

figure; hold on;
plot(mean(temp3_PO_L), 'Color', [0 0 1 1], 'LineWidth', 2); %clear temp_PO_L temp_GP_L
plot(mean(temp3_GP_L), '--', 'Color', [0 0 1 1], 'LineWidth', 2); %clear temp_PO_L temp_GP_L

plot(mean(temp3_PO_R), 'Color', [1 0 0 1], 'LineWidth', 2); %clear temp_PO_R temp_GP_R
plot(mean(temp3_GP_R), '--', 'Color', [1 0 0 1], 'LineWidth', 2); %clear temp_PO_R temp_GP_R

clear temp3_GP_R temp3_GP_L temp3_PO_R temp3_PO_L

%% tinmeg1 vs tinmeg3 BBN comparison
ylims = [0 10^-11];
xlims = [find(tlk_all_sub.tinmeg1.time == -0.300) 
         find(tlk_all_sub.tinmeg1.time ==  0.400)];

%Inhibition in top grad (!! same sensors as for tinmeg1)
topLind = find(ismember(tlk_all_sub.tinmeg3.label, topL));
topRind = find(ismember(tlk_all_sub.tinmeg3.label, topR));

%Gather tinmeg1 condition
for i = 1:numel(tlk_all_sub.tinmeg1.ID)
    temp_PO_L(i,:) = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(topLind,:);
    temp_PO_R(i,:) = tlk_all_sub.tinmeg1.PO60.PO60_90{i}(topRind,:);

    temp_GP_L(i,:) = tlk_all_sub.tinmeg1.GP60.GP60_i240{i}(topLind,:);
    temp_GP_R(i,:) = tlk_all_sub.tinmeg1.GP60.GP60_i240{i}(topRind,:);
end

%Gather tinmeg3 condition
for i = 1:numel(tlk_all_sub.tinmeg3.ID)
    temp3_PO_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.PO_00{i}(topLind,:);
    temp3_PO_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.PO_00{i}(topRind,:);

    temp3_GP_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.GPP_00{i}(topLind,:);
    temp3_GP_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.GPP_00{i}(topRind,:);
end

%RIGHT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp_PO_R), 'Color', [0 0 0 1], 'LineWidth', 1.5); %clear temp_PO_L temp_GP_L
plot(mean(temp_GP_R), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5); %clear temp_PO_R temp_GP_R

plot(mean(temp3_PO_R), 'Color', [1 0 0 1], 'LineWidth', 1.5); %clear temp_PO_R temp_GP_R
plot(mean(temp3_GP_R), '--', 'Color', [1 0 0 1], 'LineWidth', 1.5); %clear temp_PO_R temp_GP_R

ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('Non-tinnitus vs Tinnitus subjects GPIAS (Right)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

%saveas(gcf, '../analysis_output/figures/tin_vs_con_R.svg');
%close;

% LEFT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp_PO_L), 'Color', [0 0 0 1], 'LineWidth', 1.5); %clear temp_PO_L temp_GP_L
plot(mean(temp_GP_L), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5); %clear temp_PO_R temp_GP_R

plot(mean(temp3_PO_L), 'Color', [0 0 1 1], 'LineWidth', 1.5); %clear temp_PO_R temp_GP_R
plot(mean(temp3_GP_L), '--', 'Color', [0 0 1 1], 'LineWidth', 1.5); %clear temp_PO_R temp_GP_R

ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('Non-tinnitus vs Tinnitus subjects GPIAS (Left)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

%saveas(gcf, '../analysis_output/figures/tin_vs_con_L.svg');
%close;

clear temp3_GP_R temp3_GP_L temp3_PO_R temp3_PO_L temp_GP_R temp_GP_L temp_mean_R temp_mean_L temp_PO_L temp_PO_R

%% tinmeg2 vs tinmeg3 BBN

%Inhibition in top grad (!! same sensors as for tinmeg1)
topLind = find(ismember(tlk_all_sub.tinmeg3.label, topL));
topRind = find(ismember(tlk_all_sub.tinmeg3.label, topR));

%Gather tinmeg2 condition
for i = 1:numel(tlk_all_sub.tinmeg2.ID)
    temp2_PO_BBN_L(i,:) = tlk_all_sub.tinmeg2.tin0_bkg0.PO_00{i}(topLind,:);
    temp2_PO_BBN_R(i,:) = tlk_all_sub.tinmeg2.tin0_bkg0.PO_00{i}(topRind,:);

    temp2_GP_BBN_L(i,:) = tlk_all_sub.tinmeg2.tin0_bkg0.GPP_00{i}(topLind,:);
    temp2_GP_BBN_R(i,:) = tlk_all_sub.tinmeg2.tin0_bkg0.GPP_00{i}(topRind,:);
end

%Gather tinmeg 3 condition
for i = 1:numel(tlk_all_sub.tinmeg3.ID)
    temp3_PO_BBN_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.PO_00{i}(topLind,:);
    temp3_PO_BBN_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.PO_00{i}(topRind,:);

    temp3_GP_BBN_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.GPP_00{i}(topLind,:);
    temp3_GP_BBN_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg0.GPP_00{i}(topRind,:);
end

%RIGHT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp2_PO_BBN_R), 'Color', [0 0 0 1], 'LineWidth', 1.5); 
plot(mean(temp2_GP_BBN_R), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5);

plot(mean(temp3_PO_BBN_R), 'Color', [1 0 0 1], 'LineWidth', 1.5);
plot(mean(temp3_GP_BBN_R), '--', 'Color', [1 0 0 1], 'LineWidth', 1.5);

ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('BBN TinMEG2 vs TinMEG3 GPIAS (Right)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

%saveas(gcf, '../analysis_output/figures/dat2vsdat3_BBN_R.svg');
%close;

% LEFT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp2_PO_BBN_L), 'Color', [0 0 0 1], 'LineWidth', 1.5);
plot(mean(temp2_GP_BBN_L), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5);

plot(mean(temp3_PO_BBN_L), 'Color', [0 0 1 1], 'LineWidth', 1.5); 
plot(mean(temp3_GP_BBN_L), '--', 'Color', [0 0 1 1], 'LineWidth', 1.5); 
ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('BBN TinMEG2 vs TinMEG3 GPIAS (Left)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

%saveas(gcf, '../analysis_output/figures/dat2vsdat3_BBN_L.svg');
%close;

clear temp2_PO_BBN_L temp2_GP_BBN_L temp3_PO_BBN_L temp3_GP_BBN_R temp2_GP_BBN_R temp2_PO_BBN_R temp3_GP_BBN_L temp3_PO_BBN_R;

%% tinmeg2 vs tinmeg3 NBN3

%Inhibition in top grad (!! same sensors as for tinmeg1)
topLind = find(ismember(tlk_all_sub.tinmeg3.label, topL));
topRind = find(ismember(tlk_all_sub.tinmeg3.label, topR));

%Gather tinmeg2 condition
for i = 1:numel(tlk_all_sub.tinmeg2.ID)
    temp2_PO_NBN3_L(i,:) = tlk_all_sub.tinmeg2.tin0_bkg3.PO_03{i}(topLind,:);
    temp2_PO_NBN3_R(i,:) = tlk_all_sub.tinmeg2.tin0_bkg3.PO_03{i}(topRind,:);

    temp2_GP_NBN3_L(i,:) = tlk_all_sub.tinmeg2.tin0_bkg3.GPP_03{i}(topLind,:);
    temp2_GP_NBN3_R(i,:) = tlk_all_sub.tinmeg2.tin0_bkg3.GPP_03{i}(topRind,:);
end

%Gather tinmeg 3 condition
for i = 1:numel(tlk_all_sub.tinmeg3.ID)
    temp3_PO_NBN3_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg3.PO_03{i}(topLind,:);
    temp3_PO_NBN3_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg3.PO_03{i}(topRind,:);

    temp3_GP_NBN3_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg3.GPP_03{i}(topLind,:);
    temp3_GP_NBN3_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg3.GPP_03{i}(topRind,:);
end

%RIGHT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp2_PO_NBN3_R), 'Color', [0 0 0 1], 'LineWidth', 1.5); 
plot(mean(temp2_GP_NBN3_R), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5);

plot(mean(temp3_PO_NBN3_R), 'Color', [1 0 0 1], 'LineWidth', 1.5);
plot(mean(temp3_GP_NBN3_R), '--', 'Color', [1 0 0 1], 'LineWidth', 1.5);

ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('NBN3 TinMEG2 vs TinMEG3 (Right)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

saveas(gcf, '../analysis_output/figures/dat2vsdat3_NBN3_R.svg');
close;

% LEFT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp2_PO_NBN3_L), 'Color', [0 0 0 1], 'LineWidth', 1.5);
plot(mean(temp2_GP_NBN3_L), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5);

plot(mean(temp3_PO_NBN3_L), 'Color', [0 0 1 1], 'LineWidth', 1.5); 
plot(mean(temp3_GP_NBN3_L), '--', 'Color', [0 0 1 1], 'LineWidth', 1.5); 
ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('NBN3 TinMEG2 vs TinMEG3 GPIAS (Left)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

saveas(gcf, '../analysis_output/figures/dat2vsdat3_NBN3_L.svg');
close;

clear temp2_PO_NBN3_L temp2_GP_NBN3_L temp3_PO_NBN3_L temp3_GP_NBN3_R temp2_GP_NBN3_R temp2_PO_NBN3_R temp3_GP_NBN3_L temp3_PO_NBN3_R;

%% tinmeg2 vs tinmeg3 NBN8

%Inhibition in top grad (!! same sensors as for tinmeg1)
topLind = find(ismember(tlk_all_sub.tinmeg3.label, topL));
topRind = find(ismember(tlk_all_sub.tinmeg3.label, topR));

%Gather tinmeg2 condition
for i = 1:numel(tlk_all_sub.tinmeg2.ID)
    temp2_PO_NBN8_L(i,:) = tlk_all_sub.tinmeg2.tin0_bkg8.PO_08{i}(topLind,:);
    temp2_PO_NBN8_R(i,:) = tlk_all_sub.tinmeg2.tin0_bkg8.PO_08{i}(topRind,:);

    temp2_GP_NBN8_L(i,:) = tlk_all_sub.tinmeg2.tin0_bkg8.GPP_08{i}(topLind,:);
    temp2_GP_NBN8_R(i,:) = tlk_all_sub.tinmeg2.tin0_bkg8.GPP_08{i}(topRind,:);
end

%Gather tinmeg 3 condition
for i = 1:numel(tlk_all_sub.tinmeg3.ID)
    temp3_PO_NBN8_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg8.PO_08{i}(topLind,:);
    temp3_PO_NBN8_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg8.PO_08{i}(topRind,:);

    temp3_GP_NBN8_L(i,:) = tlk_all_sub.tinmeg3.tin0_bkg8.GPP_08{i}(topLind,:);
    temp3_GP_NBN8_R(i,:) = tlk_all_sub.tinmeg3.tin0_bkg8.GPP_08{i}(topRind,:);
end

%RIGHT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp2_PO_NBN8_R), 'Color', [0 0 0 1], 'LineWidth', 1.5); 
plot(mean(temp2_GP_NBN8_R), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5);

plot(mean(temp3_PO_NBN8_R), 'Color', [1 0 0 1], 'LineWidth', 1.5);
plot(mean(temp3_GP_NBN8_R), '--', 'Color', [1 0 0 1], 'LineWidth', 1.5);

ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('NBN8 TinMEG2 vs TinMEG3 (Right)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

saveas(gcf, '../analysis_output/figures/dat2vsdat3_NBN8_R.svg');
close;

% LEFT
figure('Position', [500 500 1400 500]); hold on;
plot(mean(temp2_PO_NBN8_L), 'Color', [0 0 0 1], 'LineWidth', 1.5);
plot(mean(temp2_GP_NBN8_L), '--', 'Color', [0 0 0 1], 'LineWidth', 1.5);

plot(mean(temp3_PO_NBN8_L), 'Color', [0 0 1 1], 'LineWidth', 1.5); 
plot(mean(temp3_GP_NBN8_L), '--', 'Color', [0 0 1 1], 'LineWidth', 1.5); 
ylim(ylims)
xlim(xlims)

xline(find(tlk_all_sub.tinmeg1.time == 0), '--k')

xticks([1:20:201])
xticklabels([-500:100:500])

legend({'Pulse only control', 'Gap+Pulse control', 'Pulse only tinnitus', 'Gap+Pulse tinnitus'}, 'Box','off', 'Location','northwest', 'AutoUpdate','off')

title('NBN8 TinMEG2 vs TinMEG3 GPIAS (Left)')
xlabel('Time (ms)')
ylabel('ERF Amplitude (T/cm)')

saveas(gcf, '../analysis_output/figures/dat2vsdat3_NBN8_L.svg');
close;

clear temp2_PO_NBN8_L temp2_GP_NBN8_L temp3_PO_NBN8_L temp3_GP_NBN8_R temp2_GP_NBN8_R temp2_PO_NBN8_R temp3_GP_NBN8_L temp3_PO_NBN8_R;