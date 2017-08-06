function [data, varargout] = fread(filename, varargin)
%FREAD Top-level function for reading image data from file.
%
%   DATA = FREAD(FILENAME) attempts to infer the format of the file from
%   its file extention and output the data in DATA.

[~, ~, fext] = fileparts(filename);
fext = lower(fext);

if strcmp(fext, '.tif') || strcmp(fext, '.tiff')
    [data, tags] = io.tiff.fread(filename, varargin{:});
    if nargout > 1
        varargout{1} = tags;
    end
else
    error(generatemsgid('UnknownExt'), 'Unknown file format.');
end

end
