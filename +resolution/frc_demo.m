clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'Hela_POM121_1_postprocessed.dat'));

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
nd = 1960;
% n trials
n = 20;

for i = 1:n
    


tic;
[frc_frq, frc_raw, frc_avg, frc_std, frc_num] = resolution.frccurve(coords, nd);
t = toc;
fprintf('%.2fs elapsed\n', t);

[res, frc_thr] = resolution.frcc2res(frc_frq, frc_avg);
fprintf('resolution = %.2fnm\n', res);

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
    hold on;
    plot(frc_frq, frc_thr);
    
figure('Name', 'Spurious Correlation', 'NumberTitle', 'off');
frc_num = mean(frc_num);
plot(frc_frq, frc_num);
    axis([frc_frq(1), frc_frq(end), -1, 1]);
    xlabel('Spatial Frequency (nm^{-1})');
