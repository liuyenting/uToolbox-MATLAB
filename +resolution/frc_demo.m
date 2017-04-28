clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'subarea3_frc.dat'));

fprintf(' %d samples loaded\n', size(coords, 1));

% resolution [dx, dy, dz] in nm
pxsize = [103, 103, 1000];
% magnification
mag = 5;

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
F0 = fft2(fftshift(I0p));
F1 = fft2(fftshift(I1p));
subplot(2, 2, 3);
imagesc(100*log(1+abs(fftshift(F0))));
axis image;
subplot(2, 2, 4);
imagesc(100*log(1+abs(fftshift(F1))));
axis image;

% calcluate the radial sum and their correlations

minscale = min(npx(1), npx(2));

oversmpl = 5;
frcres = linspace(0, minscale, minscale*oversmpl);
for i = 1:length(frcres)
    r = frcres(i);

    % numerator
    num = radialsum(F0 .* conj(F1), r);
    num = real(num);

    % denominator
    S1 = radialsum(abs(F0).^2, r);
    S2 = radialsum(abs(F1).^2, r);
    den = sqrt(abs(S1 .* S2));

    % result
    frcres(i) = double(num) / double(den);
    % remove NaN
    %frc_out(isnan(frc_out)) = 0;
end
figure('Name', 'FRC resolution', 'NumberTitle', 'off');
plot(frcres);
