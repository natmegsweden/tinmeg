addpath('/home/mikkel/fieldtrip/fieldtrip/')
addpath('/home/mikkel/fieldtrip/fieldtrip/external/mne')
ft_defaults



fname = '/archive/20057_working_memory/ME/NatMEG_0421/170915/memo_task-1_mc_tsss.fif';
fname = '/home/mikkel/WorkingMemory/pilot_analysis/pilot9_n-back/MEG/n_back_tica-raw.fif'


%% TEST RRIALFUN
  cfg                     = [];   
  cfg.dataset             = fname;
  cfg.trialdef.prestim    = 1;
  cfg.trialdef.poststim   = 1;
  cfg.trialdef.eventvalue = [];
  cfg.trialdef.eventtype  = 'STI101';
  cfg.trialfun            = 'ft_trialfun_neuromagSTI016fix';
  cfg                     = ft_definetrial(cfg);
  data                    = ft_preprocessing(cfg);

%% 

eve = ft_read_event(fname) %'chanindx',336)

hdr = ft_read_header(fname)

dat = ft_read_data(fname, 'chanindx', 321:338);

val = [eve.value];
smp = [eve.sample];
typ = {eve.type};
zz = find(val==0)

zztime = smp(zz)

[tf,idx] = ismember(zztime,smp)

typ = typ(zz)

sample  = unique(smp)';
latency = (sample-1)/hdr.Fs;
type    = unique(typ)';

trigarray = nan(length(sample), length(type));

for i=1:numel(sample)
  sel = find(smp==sample(i));
  for j=1:numel(sel)
      trigarray(i, strcmp(type, typ{sel(j)})) = val(sel(j));
  end
end

trigtable = array2table(trigarray, 'VariableNames', type);
trigtable = [table(sample, latency) trigtable];

writetable(trigtable, 'trigger.xls');


%% 

% Loop over files here...
dat = ft_read_data(fname, 'chanindx', 321:338); % Read the trigger data

negvals = find(dat(18,:)==-32768)
fix = dat(17,:);fix(negvals) = fix(negvals)+32768*2;

allsti = zeros(16,length(dat));

allsti = dat(1:16,:)==5;

% allsti(1,:) = find(dat(1,:)==5)
% allsti(1,:) = find(dat(2,:)==5)
% allsti(1,:) = find(dat(3,:)==5)
% allsti(1,:) = find(dat(4,:)==5)
% allsti(1,:) = find(dat(5,:)==5)
% allsti(1,:) = find(dat(6,:)==5)
% allsti(1,:) = find(dat(7,:)==5)
% allsti(1,:) = find(dat(8,:)==5)
% allsti(1,:) = find(dat(9,:)==5)
% allsti(1,:) = find(dat(10,:)==5)
% allsti(1,:) = find(dat(11,:)==5)
% allsti(1,:) = find(dat(12,:)==5)
% allsti = find(dat(13,:)==5)
% allsti = find(dat(14,:)==5)
% allsti = find(dat(15,:)==5)
% allsti = find(dat(16,:)==5)

sti101 = allsti(1,:)*2^0 + allsti(2,:)*2^1 + allsti(3,:)*2^2 + allsti(2,:)*2^3 + ...
    + allsti(5,:)*2^4 + allsti(6,:)*2^5 + allsti(3,:)*2^2 + allsti(2,:)*2^3 +

% sti101 = dat(17,:);
% sti101(sixteenON) = sti101(sixteenON)+2^15

% Can FieldTrip read this?

find(sti101==1)

% trl = [smaple_onest offset startsample code]
% [1234, 2234, -500, 1]

% cfg.filnamn = ...
% cfg.prestim = ...
% cfg.poststim = ..
% trl = ft_definetrial(cfg)
% 
% cfg.trl = trl
% dat{filenumber} = ft_preprocessing(...)

% end loop here

mydata = ft_appenddata(dat{:})
