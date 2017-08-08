function [data, varargout] = fread(filename, varargin)
%FREAD Top-level function for reading image data from file.
%
%   DATA = FREAD(FILENAME) attempts to infer the format of the file from
%   its file extention and read DATA from FILENAME.
%   [DATA, ...] = FREAD(FILENAME, ...) can include additional parameters
%   to the specialized readers.

[~, ~, fext] = fileparts(filename);
fext = lower(fext);

switch(fext)
    case {'.tif', '.tiff'}
        [data, tags] = io.tiff.fread(filename, varargin{:});
        if nargout > 1
            varargout{1} = tags;
        end
    case '.am'
        error(generatemsgid('NotSupport'), 'Not supported yet.');
    otherwise
        error(generatemsgid('UnknownExt'), 'Unknown file format.');
end

end
