clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'subarea3_frc.dat'));

% resolution [dx, dy, dz] in nm
pxsize = [103, 103, 1000];
% magnification
mag = 10;

fprintf(' %d samples loaded\n', size(coords, 1));

%% prepare the data set
fprintf('\n -- prepare the data set --\n');

% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

% estimate the output size
[npx, pxsize] = estsize(coords, pxsize, mag);

% permuted indices
permInd = randperm(size(coords, 1));
% permute the data
coords = coords(permInd, :);

tic;

I0 = resolution.binlocal(coords(2:2:end, :), npx, pxsize);
I1 = resolution.binlocal(coords(1:2:end, :), npx, pxsize);

t = toc;

fprintf(' x=%d, y=%d, z=%d\n', npx(1), npx(2), npx(3));
fprintf(' %.2fms to bin the image\n', t*1e3);

% generate Z projection
I0p = sum(I0, 3);
I1p = sum(I1, 3);

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

frc_res = resolution.frc(I0p, I1p, npx, pxsize);

figure('Name', 'FRC resolution', 'NumberTitle', 'off');

frc_frq = 0:length(frc_res)-1;
frc_frq = frc_frq / size(I0p, 1);
plot(frc_frq, frc_res);
axis([frc_frq(1), frc_frq(end), -0.5, 1]);
xlabel('Spatial Frequency (nm^{-1})');

nr = radialsum(ones(size(I0p)));
