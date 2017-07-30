function list = rmcontent(path)
%RMDIR Remove everything in the specified path.
%   Detailed explanation goes here

if isdir(path)
    content = dir(path);
    
    % ignore '.' and '..'
    nFile = numel(content);
    if nFile > 2
        content = content(3:end);
        nFile = nFile-2;
    else
        return;
    end
    
    list = cell([nFile, 1]);
    for iFile = 1:nFile
        c = content(iFile);
        
        % remove the item dependes on file or directory
        p = fullfile(path, c.name);
        if c.isdir
            rmdir(p, 's');
        else
            delete(p);
        end
        
        % save the delete path to list
        list{iFile} = p;
    end
else
    error('util:rmcontent', 'Path is not a directory.');
end

end

