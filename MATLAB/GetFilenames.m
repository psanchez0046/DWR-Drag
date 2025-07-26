function [ names_cell ] = GetFilenames(pattern, inputDataFilepath)
% Returns a cell array of filenames from a directory that match a pattern.
dinfo = dir(inputDataFilepath);             % Get directory info
all_names = {dinfo.name};                   % Extract all file names
names_cell = {};                            % Initialize output cell array
idx = 1;

for i = 1:length(all_names)
    if contains(all_names{i}, pattern)      % Use 'contains' instead of strfind
        names_cell{idx} = all_names{i};     % Add matching filename
        idx = idx + 1;
    end
end
end
