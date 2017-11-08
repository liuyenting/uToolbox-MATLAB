close all;
clearvars;

psf = 'F:\\cell_smpl.tif';

nori = 1;
nphases = 5;

%% load the data
Ipsf = io.fread(psf);
Ipsf = single(Ipsf);

%% separate the orientations and phases
Ipsf = sim.opmajor(Ipsf, nori, nphases);
volsz = size(Ipsf);
volsz = volsz(1:3);

%% z projection
Ipsf = sum(Ipsf, 3);
Ipsf = squeeze(Ipsf);

% preview
figure('Name', 'Raw Data', 'NumberTitle', 'off');
for i = 1:nphases
    subplot(1, nphases, i);
    imagesc(Ipsf(:, :, i));
        axis image;
        colormap(gray);
end

%% perform FT
Fpsf = zeros(size(Ipsf));
for i = 1:nphases
    Fpsf(:, :, i) = fftshift(fft2(ifftshift(Ipsf(:, :, i))));
end

%% calculate the power spectrum
Ppsf = zeros(size(Ipsf));
for i = 1:nphases
    Ppsf(:, :, i) = abs(Fpsf(:, :, i));
end

% preview
figure('Name', 'Power Spectrum', 'NumberTitle', 'off');
for i = 1:nphases
    subplot(1, nphases, i);
    imagesc(Ppsf(:, :, i).^0.1);
        axis image;
end