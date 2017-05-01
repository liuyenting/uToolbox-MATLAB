function compilemex(srcdir, mexname, config)
%COMPILEMEX Compile MEX file automatically.
%   Detailed explanation goes here

%% configure the compiler type
if iscudasrc(srcdir)
    compile = @mexcuda;
else
    compile = @mex;
end



end

function b = iscudasrc(srcdir)
%ISCUDASRC Verify whether the source files contain CUDA contents.

files = dir(fullfile(srcdir, '*.cu'));
b = (numel(size(files, 1)) ~= 0);
    
end
