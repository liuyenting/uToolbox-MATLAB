function imsave(I, outPath, overwrite)
%IMSAVE Save as TIFF file.
%   IMSAVE(I, OUTPATH) saves the image I to OUTPATH. If the file exists, it
%   will show error message.
%
%   I can be either a 2-D array or a 3-D array. It will save as a TIFF
%   stack if it is a 3-D array. Dimension larger then 3 is not allowed.
%
%   IMSAVE(I, OUTPATH, OVERWRITE) uses the OVERWRITE flag to indicate
%   whether target file needs to be overwritten.

if nargin == 2
    overwrite = false;
end

%% Acquire image information.
[ny, nx, nz] = size(I);
imgType = class(I);

%% Create TIFF header.
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.ImageLength = ny;
tagstruct.ImageWidth = nx;
tagstruct.RowsPerStrip = nx;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% Single pixel contains two samples if it is complex.
if isreal(I)
    tagstruct.SamplesPerPixel = 1;
else
    tagstruct.SamplesPerPixel = 2;
end

%TODO: Adjust the compression option.
tagstruct.Compression = Tiff.Compression.None;

switch imgType
    case {'uint8', 'uint16', 'uint32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'in16', 'int32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
    case {'single', 'double', 'uint64', 'int64'}
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        error('tiff:imsave:SampleFormat', 'Unknown image data type.');
end

%TODO: Identify bits per sample through imgType.
switch imgType
    case {'uint8', 'int8'}
        tagstruct.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tagstruct.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tagstruct.BitsPerSample = 32;
    case {'single'}
        tagstruct.BitsPerSample = 32;
    case {'double', 'uint64', 'int64'}
        tagstruct.BitsPerSample = 64;
    otherwise
        error('tiff:imsave:BitsPerSample', 'Unknown image data type.');
end

%% Delete existed file.
if exist(outPath, 'file')
    if overwrite
        warning('imsave:Delete', 'Remove existed output file.');
        delete(outPath);
    else
        error('tiff:imsave:Exists', 'Output file exists.');
    end
end

%% Create new TIFF file.
tiffFile = Tiff(outPath, 'a');

%% Write each layer of image to the file.
for iz = 1:nz
    tiffFile.setTag(tagstruct);
    tiffFile.write(I(:, :, iz));
    tiffFile.writeDirectory();
end
tiffFile.close();

end
