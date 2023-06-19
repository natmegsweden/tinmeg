%% Find relevant files for subject and create cell-array of file paths

% Create cell array for subjects filepaths
subpaths = cell(1);

for i = 1:numel(sub_date.ID);
    
    % Find files in subjects path with keywords specified for find_files(folder, inc_str, exc_str)
    % Note file name and path are different for different sets of data
    % colletion
    if sub_date.Exp{i} == 'tinmeg1'
        subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
        fnames = find_files(subpath, {'tinmeg1', 'tsss'}, 'ds');
    elseif sub_date.Exp{i} == 'tinmeg2'
        subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/'];
        fnames = find_files(subpath, {'tinmeg2', 'tsss'}, 'ds');
    elseif sub_date.Exp{i} == 'tinmeg3'
        subpath = [meg_data_path 'NatMEG_' char(sub_date.ID{i}) '/' char(sub_date.date{i}) '/meg/'];
        fnames = find_files(subpath, {'tinmeg3', 'tsss'}, 'ds');
    end

    for fileindex = 1:length(fnames);
        subpaths{i,1+fileindex} = [subpath char(fnames(fileindex))]; % NB! n of files differ between rows, some subjects have empty columns
    end

    subpaths{i,1} = char(sub_date.ID{i});

end

writetable(cell2table(subpaths), '../analysis_ouput/logs/included_filepaths.csv') %Write log

clear fnames subpath fileindex i


