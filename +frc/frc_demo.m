clear all; close all; %#ok<CLALL>

fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'FRC1.dat'));

% resolution [dx, dy, dz] in nm
pixelsize = [103, 103, 1000];

I = frc.binlocal(coords(:, 2:4), pixelsize, 5);
[nx, ny, nz] = size(I);

figure('Name', 'Binned', 'NumberTitle', 'off');
I_zproj = sum(I, 3);
imagesc(I_zproj);
title('Z projection');
axis image;
