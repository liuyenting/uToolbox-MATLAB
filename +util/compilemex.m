function h = compilemex(cwd, mexname, varargin)
%COMPILEMEX Compile MEX file automatically.
%   
%   CWD Current working directory

if iscompiled(cwd, mexname)
    h = str2func(mexname);
    return;
end

funcSrcPath = locatesrc(cwd, mexname);

%% configure the compiler type
if iscudasrc(funcSrcPath)
    compiler = @mexcuda;
else
    compiler = @mex;
end

%% compile the source code


end

function b = iscompiled(cwd, mexname)
%ISCOMPILED validate whether the MEX file is compiled

% append the system dependent file extension
mexname = [mexname, '.', mexext];

% search for the MEX file in path
publicPath = fullfile(cwd, mexname);
privatePath = fullfile(cwd, 'private', mexname);
b = (exist(publicPath, 'file') || exist(privatePath, 'file'));

end

function p = locatesrc(cwd, mexname)
%LOCATESRC Locate the actual source folder of target MEX function.
%
%   Note:
%   Source codes are placed under a folder with the same name as the
%   desired name for the MEX function.

p = fullfile(cwd, mexname);
if ~exist(p, 'dir')
    p = fullfile(cwd, 'private', mexname);
    if ~exist(p, 'dir')
        error('util:compilemex', 'Unable to locate the source directory.');
    end
end
    
end

function b = iscudasrc(srcdir)
%ISCUDASRC Verify whether the source files contain CUDA contents.

files = dir(fullfile(srcdir, '*.cu'));
b = (numel(size(files, 1)) ~= 0);
    
end
