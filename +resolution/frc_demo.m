clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'example_fig2a.dat'));

% resolution [dx, dy, dz] in nm
%pxsize = [103, 103, 1000];
pxsize = [103, 103, 1000];
% magnification
mag = 10;

fprintf(' %d samples loaded\n', size(coords, 1));

%% prepare the data set
fprintf('\n -- prepare the data set --\n');

% leave only XY values
coords = coords(:, 1:2);
% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

% estimate the output size
%[npx, pxsize] = estsize(coords, pxsize, mag);
%fprintf(' x=%d, y=%d, z=%d\n', npx(1), npx(2), npx(3));
%fprintf(' x=%d, y=%d\n', npx(1), npx(2));

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

% super-resolved image size
npx = [2560, 2560];
% n trials
n = 25;

[frc_raw, frc_avg, frc_std] = resolution.frccurve(coords, npx, n);

figure('Name', 'FRC resolution', 'NumberTitle', 'off');

nrs = length(frc_avg);
frc_frq = 0:nrs-1;
frc_frq = frc_frq / nrs;

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
