function extdriftcorr
%EXTDRIFTCORR Adjust the input data with drift lookup table.
%   Note: This implementation may cause JVM error after prolong execution
%   in background!

%% Retrieve file path.
[inDataFilename, inDataDir] = uigetfile('*.csv', ...
                                        'Where is the input data', ...
                                        'MultiSelect', 'off');
if ~inDataDir
    return;
else
    inDataPath = fullfile(inDataDir, inDataFilename);
end

[inCorrFilename, inCorrDir] = uigetfile('*.csv', ...
                                        'Where is the correction trajectory', ...
                                        inDataDir, ...
                                        'MultiSelect', 'off');
if ~inCorrDir
    return;
else
    inCorrPath = fullfile(inCorrDir, inCorrFilename);
end
defaultOutDataPath = createdefaultname(inDataPath);
[outDataFilename, outDataDir] = uiputfile('*.csv', ...
                                          'Save as', ...
                                          defaultOutDataPath);
if ~outDataDir
    error('local:extdriftcorr', 'Output file not specified.');
else
    outDataPath = fullfile(outDataDir, outDataFilename);
end

%% Load the lookup table and preprocess it.
tic;

% Retrieve the offset.
[~, xCol, yCol] = findlocationoffset(inCorrPath);

% Load the lookup table, ignore the header.
lut = csvread(inCorrPath, 1, 0);
lut = lut(:, [xCol, yCol]);

% Calculate the offset relative to average location.
nFrames = size(lut, 1);
lut = lut - repmat(mean(lut), nFrames, 1);

%% Load the data.
[fCol, xCol, yCol] = findlocationoffset(inDataPath);
inData = csvread(inDataPath, 1, 0);

tElapsed = toc;
fprintf('%f seconds to load the data\n', tElapsed);

%% Process the offset.
tic;

nRows = size(inData, 1);
for iRow = 1:nRows
    frame = inData(iRow, fCol);

    % Compensate the offset.
    inData(iRow, xCol) = inData(iRow, xCol) - lut(frame, 1);
    inData(iRow, yCol) = inData(iRow, yCol) - lut(frame, 2);
end

tElapsed = toc;
fprintf('%f seconds to process the data\n', tElapsed);

%% Save the result.
% Extract original header.
fid = fopen(inDataPath);
header = textscan(fid, '%s', 1, 'Delimiter', '\n');
header = char(header{1});
fclose(fid);

% Write back.
chunkSize = 10000;

tic;

fid = fopen(outDataPath, 'W');
fwrite(fid, [header, char(10)]);
buf = '';
for iRow = 1:nRows
    row = inData(iRow, :);
    buf = [buf, num2str(row(1))];
    for elem = row(2:end)
        buf = [buf, ',', num2str(elem)];
    end
    buf = [buf, char(10)];
    
    if ~mod(iRow, chunkSize)  
        fwrite(fid, buf);
        buf = '';
        fraction = iRow / nRows;
        fprintf('Saving data... (%d%% done)\n', round(100*fraction));
    end
end
fclose(fid);

tElapsed = toc;
fprintf('%f seconds to write back\n', tElapsed);

end

function outPath = createdefaultname(inPath)
%CREATEDEFAULTNAME Create default output path from input path.

[path, name, ext] = fileparts(inPath);
outPath = fullfile(path, [name, '_corr'], ext);

end

function [fCol, xCol, yCol] = findlocationoffset(fpath)
%FINDLOCATIONOFFSET Find the columns that store the X, Y location.

% Extract the header.
headers = extractcolumnheader(fpath);

% Reset the variables.
fCol = 0; xCol = 0; yCol = 0;

nHeaders = numel(headers);
for iHeaders = 1:nHeaders
    header = headers{iHeaders};
    if strncmpi(header, 'frame', 5)
        fCol = iHeaders;
    elseif strncmpi(header, 'x ', 2)
        xCol = iHeaders;
    elseif strncmpi(header, 'y ', 2)
        yCol = iHeaders;
    elseif fCol && xCol && yCol
        return;
    end
end

error('local:extdriftcorr', 'Cannot locate the X, Y column.');

end

function fields = extractcolumnheader(fpath)
%EXTRACTCOLUMNHEADER Extracts the column header from input file.

fid = fopen(fpath);

line = textscan(fid, '%s', 1, 'Delimiter', '\n');
line = char(line{1});
line = strrep(line, '"', '');
fields = strsplit(line, ',');

fclose(fid);

end
