clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'FRC1.dat'));

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

% n samples
nd = 1280;

tic;
[frcFrq, frcCrv] = resolution.frccurve(coords, nd, 'Iterations', 1);
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
