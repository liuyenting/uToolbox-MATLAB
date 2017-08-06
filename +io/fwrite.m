function fwrite(data, filename, varargin)
%FWRITE Top-level function for writing image data to file.
%
%   FWRITE(DATA, FILENAME) attempts to infer the format of the file from
%   its file extention and output DATA to FILENAME.
%   [DATA, ...] = FWRITE(FILENAME, ...) can include additional parameters
%   to the specialized writers.

[~, ~, fext] = fileparts(filename);
fext = lower(fext);

if strcmp(fext, '.tif') || strcmp(fext, '.tiff')
    io.tiff.fwrite(data, filename, varargin{:});
else
    error(generatemsgid('UnknownExt'), 'Unknown file format.');
end

end
