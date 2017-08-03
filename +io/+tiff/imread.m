function I = imread(filename, varargin)
%IMREAD
%
%   TBA
%
%   Note
%   ----
%   This function is only capable of reading simple TIFF image file.

%% parameters
p = inputParser;
addOptional(p, 'ShowWarnings', false);
parse(p, varargin{:});

isShowWarn = p.Results.ShowWarnings;

if ~isShowWarn
    % ignore warnings for unknown tags
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
end

%% pre-allocate
tiffObj = Tiff(filename, 'r');

% identify the dimension
nCol = tiffObj.getTag('ImageWidth');
nRow = tiffObj.getTag('ImageLength');

%TODO allow multiple samples per pixel
nDepth = tiffObj.getTag('SamplesPerPixel');
assert(nDepth == 1);

nLayer = 0;
while true
    nLayer = nLayer+1;
    if tiffObj.lastDirectory()
        break;
    else
        tiffObj.nextDirectory();
    end
end
% index values are one-based
tiffObj.setDirectory(1);

% identify the data type
bps = tiffObj.getTag('BitsPerSample');
if bps == 8
    dataType = 'uint8';
elseif bps == 16
    dataType = 'uint16';
elseif bps == 32
    dataType = 'single';
else
    error(generatemsgid('UnknownType'), 'Unknown pixel data type.');
end

% create the array
I = zeros([nRow, nCol, nLayer, nDepth], dataType);

%% load the data
nLayer = 0;
while true
    nLayer = nLayer+1;
    if tiffObj.lastDirectory()
        break;
    else
        tiffObj.nextDirectory();
    end
end
% index values are one-based
for i = 1:nLayer
    tiffObj.setDirectory(i);
    I(:, :, i, :) = tiffObj.read();
end

%% clean up 
% release the file
tiffObj.close();

if ~isShowWarn
    % re-enable the warning
    warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
end

end