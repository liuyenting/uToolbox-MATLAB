close all;
clearvars -except data;

import matlab.*;

% filePath = fullfile(userpath, 'frc_test_data', 'usaf1951', 'usaf1951_cam100nm_dp5um_fit.csv');
filePath = 'F:\lu\cell6_frc\90.csv';

% start the diary
consolelogger('start', matlab.chfext(filePath, 'txt'));

%% loading the data
fprintf('\n -- loading data --\n');

forceReload = true;
fprintf('path = "%s"\n', filePath);

tic;

% determine the loader
[~, ~, fext] = fileparts(filePath);
if strcmp(fext, '.csv')
    fprintf('loading CSV file\n');

    % load the header
    fid = fopen(filePath);
    header = fgetl(fid);
    % split by comma
    header = strsplit(header, ',');
    % remove quotes
    header = strrep(header, '"', '');
    fclose(fid);

    fprintf('%d columns in the dataset\n', length(header));

    if exist('data', 'var') && ~forceReload
        warning('resolution:frc_demo', 'Using preloaded data.');
    else
        % load the data
        data = csvread(filePath, 1, 0);
    end
    fprintf('... %d samples loaded\n', size(data, 1));
elseif strcmp(fext, '.dat')
    fprintf('loading DAT file\n');
    coords = dlmread(filePath);
    fprintf('... %d samples loaded\n', size(coords, 1));
else 
    error('resolution:frc_demo', 'Unknown input file type.');
end

t = toc;
fprintf('%.2fs elapsed\n', t);

%% prepare the data set
fprintf('\n -- prepare the data set --\n');

if strcmp(fext, '.csv')
    % find the indices
    xyIndex = findcol(header, {'x', 'y'});
    uncertaintyIndex = findcol(header, {'uncertainty_xy'});
    if isempty(xyIndex)
        error('resolution:frc_demo', ...
              'Unable to locate the coordinate columns.');
    end
    if isempty(uncertaintyIndex)
        error('resolution:frc_demo', ...
              'Unable to locate radial uncertainty column.');
    end

    % extract the data
    coords = data(:, xyIndex);
    uncertainty = data(:, uncertaintyIndex);
elseif strcmp(fext, '.dat')
    coords = coords(:, 2:3);
    uncertainty = [];
end

% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

% super-resolved resolution [nm]
res = 10;

tic;
% [frcFrq, frcCrv, frcSpu] = resolution.frccurve(coords, res, uncertainty, ...
%                                                'Iterations', 20);
[frcFrq, frcCrv] = resolution.frccurve(coords, res, 'Iterations', 20);                                           
t = toc;
fprintf('%.2fs elapsed\n', t);

hFrcCrv = figure('Name', 'FRC resolution', 'NumberTitle', 'off');
plot(frcFrq, frcCrv);
    axis tight;
    xlim([frcFrq(1), frcFrq(end)]);
    yl = ylim; yl(2) = 1; ylim(yl); % force the max scale to 1
    xlabel('Spatial Frequency (nm^{-1})');
    ylabel('FRC');

%% convert to resolution
fprintf('\n -- convert to resolution --\n');

[res, frcThr] = resolution.frcc2res(frcFrq, frcCrv);
if isinf(res)
    fprintf('unable to solve the resolution\n');
    return;
else
    fprintf('resolution = %.2fnm\n', res);
    
    figure(hFrcCrv);
    hold on;
    plot(frcFrq, frcThr);
    hold off;
end

% hFrcSpu = figure('Name', 'Spurious Correlation', 'NumberTitle', 'off');
% plot(frcFrq, frcSpu);
%     axis tight;
%     xlim([frcFrq(1), frcFrq(end)]);
%     xlabel('Spatial Frequency (nm^{-1})');
%     ylabel('log_{10}FRC numerator');
    
consolelogger('stop');
