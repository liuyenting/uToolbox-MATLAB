clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'subarea3_frc.dat'));

fprintf(' %d samples loaded\n', size(coords, 1));

% resolution [dx, dy, dz] in nm
pxsize = [103, 103, 1000];
%pxsize = [0.89, 0.81, 10000];
% magnification
mag = 10;

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
F0 = fft2(fftshift(I0p));
F1 = fft2(fftshift(I1p));
subplot(2, 2, 3);
imagesc(100*log(1+abs(fftshift(F0))));
axis image;
subplot(2, 2, 4);
imagesc(100*log(1+abs(fftshift(F1))));
axis image;

drawnow;

% numerator
frc_num = radialsum(F0 .* conj(F1));
frc_num = real(frc_num);

% denominator
S1 = radialsum(abs(F0).^2);
S2 = radialsum(abs(F1).^2);
frc_den = sqrt(abs(S1 .* S2));

% result
frc_res = double(frc_num) ./ double(frc_den);
% remove NaN
frc_res(isnan(frc_res)) = 0;

figure('Name', 'FRC resolution', 'NumberTitle', 'off');

frc_frq = 0:length(frc_res)-1;
frc_frq = frc_frq / size(I0p, 1);
plot(frc_frq, frc_res);
axis([frc_frq(1), frc_frq(end), -0.5, 1]);
xlabel('Spatial Frequency (nm^{-1})');

nr = radialsum(ones(size(I0p)));
