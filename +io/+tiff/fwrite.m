function fwrite(data, filename, varargin)
%FWRITE Save as TIFF file.
%
%   FWRITE(DATA, FILENAME) saves the image DATA to FILENAME. If the file 
%   exists, it will show error message.
%
%   I can be either a 2-D array or a 3-D array. It will save as a TIFF
%   stack if it is a 3-D array. Dimension larger then 3 is not allowed.
%
%   FWRITE(I, OUTPATH, OVERWRITE) uses the OVERWRITE flag to indicate
%   whether target file needs to be overwritten.
%
%   Note
%   ----
%   This function is aimed for microscopy image types, which are grayscale
%   by nature with mixing of IEEE floating point image formats.

%% parameters
p = inputParser;
addParameter(p, 'Overwrite', false);
addParameter(p, 'Compression', 'none', @validCompression);
parse(p, varargin{:});

overwrite = p.Results.Overwrite;
compression = p.Results.Compression;

% image dimension
[ny, nx, nz] = size(data);
dataType = class(data);

%% create TIFF header
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.ImageLength = ny;
tags.ImageWidth = nx;
tags.RowsPerStrip = nx;
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% complex data has two samples per pixel
if isreal(data)
    tags.SamplesPerPixel = 1;
else
    tags.SamplesPerPixel = 2;
end

switch compression
    case 'none'
        tags.Compression = Tiff.Compression.None;
    case 'lzw'
        tags.Compression = Tiff.Compression.LZW;
end

switch dataType
    case {'uint8', 'uint16', 'uint32'}
        tags.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'in16', 'int32'}
        tags.SampleFormat = Tiff.SampleFormat.Int;
    case {'single', 'double', 'uint64', 'int64'}
        tags.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        error(generatemsgid('UnknownDateType'), ...
              'Unknown image data type.');
end

switch dataType
    case {'uint8', 'int8'}
        tags.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tags.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tags.BitsPerSample = 32;
    case {'single'}
        tags.BitsPerSample = 32;
    case {'double', 'uint64', 'int64'}
        tags.BitsPerSample = 64;
    otherwise
        error(generatemsgid('UnknownDateType'), ...
              'Unknown image data type.');
end

% addtional info
tags.Software = 'MATLAB uToolbox';

%% Delete existed file.
if exist(filename, 'file')
    if overwrite
        warning(generatemsgid('Delete'), 'Remove existed output file.');
        delete(filename);
    else
        error(generatemsgid('FileExists'), 'Output file exists.');
    end
end

%% create new file
tiffObj = Tiff(filename, 'a');

%% write each layer of image to the file.
for iz = 1:nz
    tiffObj.setTag(tags);
    if tags.SamplesPerPixel == 2
        re = real(data(:, :, iz));
        im = imag(data(:, :, iz));
        % Real and imaginary data are separated to different channels.
        tiffObj.write(cat(4, re, im));
    else
        tiffObj.write(data(:, :, iz));
    end
    tiffObj.writeDirectory();
end
tiffObj.close();

end

function isValid = validCompression(scheme)
%VALIDCOMPRESSION Validates whether the compression scheme is supported.
%
%   ISVALID = VALIDCOMPRESSION(SCHEME) compares SCHEME to an internal list.
%
%   Note
%   ----
%   Though tifflib from MATLAB uses libtiff, which in turns support all the
%   defined compression alogorithm, for the sake of simplicity, only LZW is
%   supported for now.

switch(scheme)
    case {'none', 'lzw'}
        isValid = true;
    otherwise
        isValid = false;
end

end
