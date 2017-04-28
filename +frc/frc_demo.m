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
[npx, pxsize] = estsize(coords, pxsize, mag)

% permuted indices
permInd = randperm(size(coords, 1));
% permute the data
coords = coords(permInd, :);

tic;

I0 = frc.binlocal(coords(2:2:end, :), npx, pxsize);
I1 = frc.binlocal(coords(1:2:end, :), npx, pxsize);

t = toc;

fprintf(' x=%d, y=%d, z=% d\n', npx(1), npx(2), npx(3));
fprintf(' %.2fms to bin the image\n', t*1e3);

figure('Name', 'Binned', 'NumberTitle', 'off');
I0p = sum(I0, 3);
I1p = sum(I1, 3);
subplot(2, 2, 1);
imagesc(I0p);
axis image;
subplot(2, 2, 2);
imagesc(I1p);
axis image;

%% masking spatial domain
fprintf('\n -- masking spatial domain --\n');

% generate Tukey window
mask = tukeywin2(npx, 8);

% mask the binned images
I0p = I0p .* mask;
I1p = I1p .* mask;

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

% acquire the FFT
F0p = fft2(fftshift(I0p));
F1p = fft2(fftshift(I1p));
subplot(2, 2, 3);
imagesc(100*log(1+abs(fftshift(F0p))));
axis image;
subplot(2, 2, 4);
imagesc(100*log(1+abs(fftshift(F1p))));
axis image;

% calcluate the radial sum and their correlations
