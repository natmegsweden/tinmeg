%% Specify conditions and event triggers

%all conditions
conditions = ({'tin0_bkg0', 'tin0_bkg3', 'tin0_bkg8', ...
               'tin3_bkg0', 'tin3_bkg3', 'tin3_bkg8', ...
               'tin8_bkg0', 'tin8_bkg3', 'tin8_bkg8'});

%Save structure of conditions
%save('../mat_data/conditions.mat', 'conditions');

%Structure for triggers and labels
cond = struct();

%trigger at pulse onset
cond.tin0_bkg0trig   = [1 2 4 8 16 32];
cond.tin0_bkg0label  = ({'GPP_00', 'GPG_00', 'PO_00', 'GO_00', 'PPP_00', 'PPG_00'});
cond.tin0_bkg3trig   = [49 50 52 56];
cond.tin0_bkg3label  = ({'GPP_03', 'GPG_03', 'PO_03', 'GO_03'});
cond.tin0_bkg8trig   = [33 34 36 40];
cond.tin0_bkg8label  = ({'GPP_08', 'GPG_08', 'PO_08', 'GO_08'});

cond.tin3_bkg0trig   = [193 194 196 200];
cond.tin3_bkg0label  = ({'GPP_30', 'GPG_30', 'PO_30', 'GO_30'});
cond.tin3_bkg3trig   = [241 242 244 248];
cond.tin3_bkg3label  = ({'GPP_33', 'GPG_33', 'PO_33', 'GO_33'});
cond.tin3_bkg8trig   = [225 226 228 232];
cond.tin3_bkg8label  = ({'GPP_38', 'GPG_38', 'PO_38', 'GO_38'});

cond.tin8_bkg0trig   = [129 130 132 136];
cond.tin8_bkg0label  = ({'GPP_80', 'GPG_80', 'PO_80', 'GO_80'});
cond.tin8_bkg3trig   = [177 178 180 184];
cond.tin8_bkg3label  = ({'GPP_83', 'GPG_83', 'PO_83', 'GO_83'});
cond.tin8_bkg8trig   = [161 162 164 168];
cond.tin8_bkg8label  = ({'GPP_88', 'GPG_88', 'PO_88', 'GO_88'});

%save('../mat_data/cond_tinmeg2.mat', 'cond');

%% Read subject list

%Readtable of subjects (as string)
sub_date = readtable('../sub_date_tinmeg2.txt', 'Format', '%s%s');

disp(['Number of subjects in table is ' num2str(height(sub_date))])