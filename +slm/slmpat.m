function P = slmpat
%SLMPAT Generates the SLM template.

clearvars; close all;

%% Parameters.

tic;

%
% SLM
%
% SLM dimension, [width, height] (px).
slmDim = [2048, 2048];
% Pixel size (um).
pxSize = 8;

%
% Optics
%
% Wavelength (nm).
wavelength = 488;
% Magnification.
mag = 116.9;
% Aperture diameter (mm).
naID = 0.325; % ?
naOD = 0.4;

%
% Pattern
%
% Oversampling ratio.
ovrSmpl = 1;
% Number of repetitive patterns.
nBeam = 10;
% Spacing (um).
spacing = 1.35;
% Tilt (rad).
tilt = 0.01;
% Clamp threshold, [0, 1].
clamp = 0.15;

%% Calculate the spatial and frequency resolution.
nSample = slmDim(1) * ovrSmpl;
ds = pxSize / (mag * ovrSmpl); % TODO: error?
df = (2*pi / ds) / nSample;

%% Generate the grids (um).
hSlmDim = slmDim/2;
% Generate the range vectors for calculation.
cvx = linspace(-hSlmDim(1), hSlmDim(1), nSample);
cvy = linspace(-hSlmDim(2), hSlmDim(2), nSample);

svx = cvx * ds;
svy = cvy * ds;
fvx = cvx * df;
fvy = cvy * df;

clear cvx cvy;

% Generate the grids using the range vectors.
% Spatial grid.
[sx, sy] = meshgrid(svx, svy);
% Frequency grid.
[fx, fy] = meshgrid(fvx, fvy);

clear svx svy fvx fvy;

% Actual SLM pixel grid.
rvx = linspace(-hSlmDim(1), hSlmDim(1), slmDim(1)) * ds;
rvy = linspace(-hSlmDim(2), hSlmDim(2), slmDim(2)) * ds;
[rx, ry] = meshgrid(rvx, rvy);

clear rvx rvy;

%% Generate the mask;
% NA radius (um).
naOR = naOD; % TODO: fix this, understant the def.
naIR = naID;
% Wavelength (um).
wavelength = wavelength / 1e3;

% Distance vector from center of the SLM in frequency domain.
dist = sqrt(fx.^2 + fy.^2);
% Frequency span of the aperture.
naOFS = naOR * (2*pi / wavelength);
naIFS = naIR * (2*pi / wavelength);
% Region inside the ring is zero.
mask = (dist < naOFS) & (dist > naIFS);

%% Populate the ideal bessel pattern in frequency space.
% Rotated frequency grid.
rotFreqGrid = fx*cos(tilt) + fy*sin(tilt);

% Frequency components should be confined in an annular region.
fPat = mask;
% Start calculating the multi-beam pattern.
for iBeam = 1:nBeam-1
    interval = spacing * iBeam;
    newField = exp(1i * rotFreqGrid*interval) + ...
               exp(-1i * rotFreqGrid*interval);
    % Mask the new field and superimpose it.
    fPat = fPat + mask.*newField;
end
% % Apply spatial shift (in frequency domain).
% shFreqGrid = fx*shiftX + fy*shiftY;
% field = field .* exp(1i * shFreqGrid);

%% Calculate the ideal result after a lens.
% Pass through a lens, therefore, an FT.
sPat = fftshift(fft2(ifftshift(fPat)));

% Retrieve the electric field and resize to [-1, 1].
eField = real(sPat);
eField = eField / max(eField(:));

% Calculate the intensity field.
iField = abs(eField).^2;

figure('Name', 'Ideal Patterns', 'NumberTitle', 'off');
subplot(1, 2, 1); 
    imagesc(fPat);
    axis equal off;
    title('Prior to the Mask');
subplot(1, 2, 2); 
    imagesc(iField); 
    axis equal off;
    title('Posterior to the Lens');
    
%% Binarize the result.
% Downsample the data.
if ovrSmpl > 1
    slmPat = interp2(sx, sy, eField, rx, ry);
else
    slmPat = eField;
end
% Clamp the value.
slmPat(abs(slmPat) < clamp) = 0;
% Show the negative region only.
slmPat = slmPat < 0;

figure('Name', 'SLM Pattern', 'NumberTitle', 'off');
    imagesc(slmPat);

toc;

%% Estimate the result.

P = slmPat;

end

