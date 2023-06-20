%% Find sensor with highest amplitude ERF in timlockeds

%Load if not already in workspace
if exist('tlk_all_sub','var') == 0
tlk_all_sub = load('../processed_data/timelockeds/aggregated/tlk_all_sub_cmb.mat');
tlk_all_sub = tlk_all_sub.tlk_all_sub;
end

%Load sensor information
load('../input_vars/sensors.mat');

% Find gradiometers with highest response in PO60_90 TOI
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