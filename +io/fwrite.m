function fwrite(data, filename, varargin)
%FWRITE Top-level function for writing image data to file.
%
%   TBA

[~, ~, fext] = fileparts(filename);
fext = lower(fext);

if strcmp(fext, '.tif') || strcmp(fext, '.tiff')
    io.tiff.fwrite(data, filename, varargin);
else
    error(generatemsgid('UnknownExt'), 'Unknown file format.');
end

end
