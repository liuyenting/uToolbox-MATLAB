function fwrite(data, filename, varargin)
%FWRITE Top-level function for writing image data to file.
%
%   FWRITE(DATA, FILENAME) attempts to infer the format of the file from
%   its file extention and output DATA to FILENAME.
%   [DATA, ...] = FWRITE(FILENAME, ...) can include additional parameters
%   to the specialized writers.

[~, ~, fext] = fileparts(filename);
fext = lower(fext);

switch(fext)
    case {'.tif', '.tiff'}
        io.tiff.fwrite(data, filename, varargin{:});
    case '.am'
        error(generatemsgid('NotSupport'), 'Not supported yet.');
    otherwise
        error(generatemsgid('UnknownExt'), 'Unknown file format.');
end

end
