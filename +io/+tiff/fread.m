function [data, varargout] = fread(filename, varargin)
%FREAD Read a TIFF file.
%
%   TBA
%
%   Note
%   ----
%   This function is aimed for microscopy image types, which are grayscale
%   by nature with mixing of IEEE floating point image formats.

%% parameters
p = inputParser;
addOptional(p, 'Warnings', false);
parse(p, varargin{:});

showWarnings = p.Results.Warnings;

if ~showWarnings
    % ignore warnings for unknown tags
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
end

%% read the tags
% 4th parameters set to 0 to read all the IDF, in order to determine how 
% many layers in this file.
tags = matlab.io.internal.imagesci.tifftagsread(filename, 0, 0, 0);
nLayer = numel(tags);
% use the info from first layer ONLY
tags = tags(1);

%% parse the format
tiffObj = Tiff(filename, 'r');

% dimension
nCol = tags.Width;
nRow = tags.Height;

% pixel format
nDepth = tags.SamplesPerPixel;
assert(nDepth == 1);

% data type
switch(tags.BitsPerSample)
    case 8
        dataType = 'uint8';
    case 16
        dataType = 'uint16';
    case 32
        dataType = 'single';
    otherwise
        error(generatemsgid('UnknownType'), 'Unknown pixel data type.');
end

%% pre-allocation
% create the array
data = zeros([nRow, nCol, nLayer, nDepth], dataType);

%% load the data
% index values are one-based
for i = 1:nLayer
    tiffObj.setDirectory(i);
    data(:, :, i, :) = tiffObj.read();
end

%% clean up 
% release the file
tiffObj.close();

if ~showWarnings
    % re-enable the warning
    warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
end

%% output the result
if nargout > 1
    varargout{1} = tags;
end
% additional output arguments are ignored

end