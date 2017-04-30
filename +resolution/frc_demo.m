clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'subarea3_frc.dat'));

% resolution [dx, dy, dz] in nm
pxsize = [103, 103];

fprintf('%d samples loaded\n', size(coords, 1));

%% prepare the data set
fprintf('\n -- prepare the data set --\n');

% leave only XY values
coords = coords(:, 1:2);
% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

% super-resolved image size
npx = [1960, 1960];
% n trials
n = 25;

tic;
[frc_frq, frc_raw, frc_avg, frc_std] = resolution.frccurve(coords, npx, n);
t = toc;
fprintf('%.2fs elapsed\n', t);

figure('Name', 'FRC resolution', 'NumberTitle', 'off');

subplot(2, 1, 1);
    plot(frc_frq, frc_raw);
        axis([frc_frq(1), frc_frq(end), -1, 1]);
        xlabel('Spatial Frequency (nm^{-1})');
        title('Raw');

subplot(2, 1, 2);
    plot(frc_frq, frc_avg);
        axis([frc_frq(1), frc_frq(end), -1, 1]);
        xlabel('Spatial Frequency (nm^{-1})');
        title('Averaged');
    hold on;
    errorbar(frc_frq, frc_avg, frc_std);
