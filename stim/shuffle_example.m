
% trigtab = readtable('AudioFile stimuli and playlist\trigtable.xlsx');
% 
% %Cell of 27 filenames
% fnames = table2cell(trigtab(2:28,1));

%Filenames of the 27 different conditions
% tin_0 = no "tinnitus", tin3 = simulated tinnitus at 3kHz, tin8 = * at 8kHz
% bkg0 = BBN background, bkg3 = NBN background CF 3kHz, bkg8 = *CF 8kHz
% _PO = Pulse only, GO = Gap only, GP = Gap + Pulse

fnames = {'tin0_bkg0_GO';
    'tin0_bkg0_PO';
    'tin0_bkg0_GP';
    'tin0_bkg3_GO';
    'tin0_bkg3_PO';
    'tin0_bkg3_GP';
    'tin0_bkg8_GO';
    'tin0_bkg8_PO';
    'tin0_bkg8_GP';
    'tin3_bkg0_GO';
    'tin3_bkg0_PO';
    'tin3_bkg0_GP';
    'tin3_bkg3_GO';
    'tin3_bkg3_PO';
    'tin3_bkg3_GP';
    'tin3_bkg8_GO';
    'tin3_bkg8_PO';
    'tin3_bkg8_GP';
    'tin8_bkg0_GO';
    'tin8_bkg0_PO';
    'tin8_bkg0_GP';
    'tin8_bkg3_GO';
    'tin8_bkg3_PO';
    'tin8_bkg3_GP';
    'tin8_bkg8_GO';
    'tin8_bkg8_PO';
    'tin8_bkg8_GP'};

%%

%repetitions of conditions within block
reps = 5;

%Create two vectors [1 2 3] and shuffle order
randvec = 1:3;
randvec = randvec(randperm(3));

randvec2 = 1:3;
randvec2 = randvec2(randperm(3));

%For each of three "tinnitus" conditions (i.e, no, low or high)
for i = 1:3
    
    %Assign a block of 9 rows to a cell
    tinidx{1} = [1:9];
    tinidx{2} = [10:18];
    tinidx{3} = [19:27];
    
    %Select 9 rows based on the shuffled order in randvec
    if i == 1
        tinblock = fnames(tinidx{randvec(1)});
    elseif i == 2
        tinblock = fnames(tinidx{randvec(2)});
    elseif i == 3
        tinblock = fnames(tinidx{randvec(3)});
    end
        
        %For each of three background noises (i.e. BBN, NBN_3, NBN_8)
        for ii = 1:3
            
            %Assign a block of 3 rows to a cell
            bkgidx{1} = [1:3];
            bkgidx{2} = [4:6];
            bkgidx{3} = [7:9];

            %Select 3 rows based on shuffled order of randvec2
            if ii == 1
                bkgblock = tinblock(bkgidx{randvec2(1)});
            elseif ii == 2
                bkgblock = tinblock(bkgidx{randvec2(2)});
            elseif ii == 3
                bkgblock = tinblock(bkgidx{randvec2(3)});
            end
            
            %Copy trials to number of repetitions sepcified
            bkgblock = repmat(bkgblock, reps, 1);
            
            %Create shuffled order of index for trials
            stimorder = 1:numel(bkgblock);
            stimorder = stimorder(randperm(length(stimorder)));
                
                %List trials according to shuffled order
                for iii = 1:numel(bkgblock)
                    disp(bkgblock{stimorder(iii)})
                end
    end
end
        
        
        