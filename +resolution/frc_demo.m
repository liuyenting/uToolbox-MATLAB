close all;
clearvars -except data;

%% loading the data
fprintf('\n -- loading data --\n');


filePath = fullfile(userpath, '0421Neuron1_postprocessed.csv');
fprintf('path = "%s"\n', filePath);

tic;

% load the header
fid = fopen(filePath);
header = fgetl(fid);
% split by comma
header = strsplit(header, ',');
% remove quotes
header = strrep(header, '"', '');
fclose(fid);

fprintf('%d columns in the dataset\n', length(header));

if exist('data', 'var')
    warning('resolution:frc_demo', ...
            'Using preloaded data.');
else
    % load the data
    data = csvread(filePath, 1, 0);
end
fprintf('... %d samples loaded\n', size(data, 1));

t = toc;
fprintf('%.2fs elapsed\n', t);


% resolution [dx, dy, dz] in nm
pxsize = [103, 103];

%% prepare the data set
fprintf('\n -- prepare the data set --\n');

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
% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

% n samples
nd = 1280;

tic;
[frcFrq, frcCrv] = resolution.frccurve(coords, nd, uncertainty, ...
                                       'Iterations', 5);
t = toc;
fprintf('%.2fs elapsed\n', t);

%[res, frc_thr] = resolution.frcc2res(frc_frq, frc_avg);
%fprintf('resolution = %.2fnm\n', res);

figure('Name', 'FRC resolution', 'NumberTitle', 'off');

plot(frcFrq, frcCrv);
    axis([frcFrq(1), frcFrq(end), -0.5, 1]);
    xlabel('Spatial Frequency (nm^{-1})');
hold on;
%plot(frcFrq, frcThr);
