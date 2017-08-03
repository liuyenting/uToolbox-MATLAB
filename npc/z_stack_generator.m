clear all

% The directory that store the TIFFs.
target = 'G:\yeast live image\WT_mtGFP_live04_5_R3D_D3D_TIFFS - Copy\';

list = dir(target);
file_counts = size(list);

% Filter out unique names.
names = cell(file_counts);
name_counts = 1;
for n = 2:file_counts;
    % Skip the directories.
    if list(n).isdir
        continue;
    end

    % Pull out the file name only.
    [~, name, ~] = fileparts(list(n).name);
    
    % Remove the Z marker.
    new_name = strsplit(name, '_z');
    new_name = new_name(1);

    % Don't store same names to the list, then no need for sorting.
    if name_counts > 1
        tmp_name = names{name_counts-1};
        if strcmp(new_name, tmp_name)
            continue;
        end
    end
    
    names{name_counts} = new_name;
    name_counts = name_counts+1;
end

% Remove empty cells.
names = names(~cellfun('isempty', names));

% Initiate the waitbar.
[total_files, ~] = size(names);
current = 0;
waitbar(current/total_files);

% Cycle through the list to look for files.
for template = names'
    new_file = strcat(target, template{:}, '.tif');
    new_file = new_file{:};
    
    % Get the associated files.
    tmp_file = strcat(target, template{:}, '*.tif');
    tmp_file = tmp_file{:};
    images = dir(tmp_file);
    
    % Create new file using the template name.
    first_time = 1;
    for image = images'
        tmp_file = strcat(target, image.name);
        
        % Read the file.
        dummy = imread(tmp_file);
        % NOTE: Maybe pre-write is not a requirement.
        if first_time
            % Create new file for the first image.
            first_time = 0;
            imwrite(dummy, new_file);
        else
            imwrite(dummy, new_file, 'writemode', 'append');
        end
    end
    
    current = current+1;
    waitbar(current/total_files);
end