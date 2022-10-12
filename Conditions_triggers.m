%% Specify conditions and event triggers

%all conditions
conditions = ({'PO60', 'PO70', 'GP60', 'GP70', 'GO'});
%Save structure of conditions
save('../mat_data/conditions.mat', 'conditions');

%Structure for triggers and labels
cond = struct();

%trigger at pulse onset
cond.PO60trig   = [40968 36872 34824 33800 33288 33032];
cond.PO60label  = ({'PO60_70', 'PO60_75', 'PO60_80', 'PO60_85', 'PO60_90', 'PO60_95'});
cond.PO70trig   = [36876 34828 33804 33292 33036];
cond.PO70label  = ({'PO70_75', 'PO70_80', 'PO70_85', 'PO70_90', 'PO70_95'});
cond.GP60trig   = [49800 49736 49704 49688];
cond.GP60label  = ({'GP60_i0', 'GP60_i60', 'GP60_i120', 'GP60_i240'});
cond.GP70trig   = [49804 49740 49708 49692];
cond.GP70label  = ({'GP70_i0', 'GP70_i60', 'GP70_i120', 'GP70_i240'});

%trigger at gap onset (gap only)
cond.GOtrig     = [16386 16390];
cond.GOlabel    = ({'GO_60', 'GO_70'});

save('../mat_data/cond.mat', 'cond');

%ISI pulse trigger (gaptrig = pulsetrig - 6)
% A_C60_i0: 49800
% A_C60_i60: 49736
% A_C60_i120: 49704
% A_C60_i240: 49688

% A_70_i0: 49804
% A_70_i60: 49740
% A_70_i120: 49708
% A_70_i240: 49692

%% Read subject list

%Readtable of subjects (as string)
sub_date = readtable('../sub_date.txt', 'Format', '%s%s');

disp(['Number of subjects in table is ' num2str(height(sub_date))])