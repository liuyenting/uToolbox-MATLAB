function poolsize = startparallel(varargin)
%STARTPARALLEL Launch worker pool if capable.
%
%   TBA

if ~util.license.islicensed('Parallel Computing Toolbox')
    warning('util:startparallel', ...
            'Parallel Computing Toolbox is required! Ignore request.');
    return;
end

p = inputParser;
addParameter(p, 'Profile', parallel.defaultClusterProfile);
addParameter(p, 'NumWorkers', feature('numCores'), ...
             @(x) ((rem(x, 1) == 0) && (x > 0)));
parse(p, varargin{:});

profile = p.Results.Profile;
nw = p.Results.NumWorkers;

% do not create a new pool when probing
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolsize = 0;
else
    restart = true;
    % verify the profile and the size match the requirement
    if ~strcmp(util.parallel.currprofile, profile)
        warning('util:parallel:startparallel', 'Profile mismatch.');
    else    
        poolsize = poolobj.NumWorkers;
        if poolsize ~= nw
            warning('util:parallel:startparallel', 'Pool size mismatch.');
        else
            restart = false;
        end
    end
    
    if restart
        delete(poolobj);
        poolsize = 0;
    end
end

if poolsize == 0
    poolobj = parpool(profile, nw);
end
poolsize = poolobj.NumWorkers;

end
