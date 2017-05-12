close all;
clearvars -except data;

%% loading the data
fprintf('\n -- loading data --\n');


filePath = fullfile(userpath, '0605cell3_driftcorrected.csv');
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

% coords = dlmread(fullfile(userpath, 'Hela_POM121_1_postprocessed.dat'));
% coords = coords(:, 1:2);
% uncertainty = [];

% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

% n samples
nd = 2560;

tic;
% [frcFrq, frcCrv, frcSpu] = resolution.frccurve(coords, nd, uncertainty, ...
%                                                'Iterations', 5);
[frcFrq, frcCrv] = resolution.frccurve(coords, nd, 'Iterations', 10);                                           
t = toc;
fprintf('%.2fs elapsed\n', t);

[res, frcThr] = resolution.frcc2res(frcFrq, frcCrv);
fprintf('resolution = %.2fnm\n', res);

figure('Name', 'FRC resolution', 'NumberTitle', 'off');
% subplot(2, 1, 1);
    plot(frcFrq, frcCrv);
        axis tight;
        xlim([frcFrq(1), frcFrq(end)]);
        yl = ylim; yl(2) = 1; ylim(yl); % force the max scale to 1
        xlabel('Spatial Frequency (nm^{-1})');
        ylabel('FRC');
    hold on;
    plot(frcFrq, frcThr);
% subplot(2, 1, 2);
%     plot(frcFrq, frcSpu);
%         axis tight;
%         xlim([frcFrq(1), frcFrq(end)]);
%         xlabel('Spatial Frequency (nm^{-1})');
%         ylabel('log_{10}FRC numerator');
