clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'FRC1.dat'));

% resolution [dx, dy, dz] in nm
pxsize = [103, 103, 1000];
% magnification
mag = 5;

%% prepare the data set
% offset back to the origin and drop the t-axis
coords = offsetorigin(coords(:, 2:4));

% estimate the output size
[npx, pxsize] = estsize(coords, pxsize, mag);

mask = tukeywin2(npx, 8);
figure();
imagesc(mask);
colormap(gray);
