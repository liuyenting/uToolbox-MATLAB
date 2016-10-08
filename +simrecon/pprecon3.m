function pprecon3
%PPRECON3 Reconstruct the 3D SIM dataset using plane-by-plane method.

clearvars; close all;

MSG_ID = 'simrecon:pprecon3';

%% Preset the parameters.
inDir = '~/Documents/SIM/RAWcell1';
psfFilename = '~/Documents/SIM/ygbead_zp1um_NAp55nap44_ExpPsf.tif';

[path, name, ~] = fileparts(inDir);
outDir = fullfile(path, [name, '_recon']);

% === Data Sets ===
% Initial position (um).
posInit = 100;
% Position increment (nm).
posDelta = 250;
% Maximum of orientation types, IT STARTS FROM 0.
oriMax = 1;
% Maximum of illumination types.
illMax = 5;

% === Optics ===
% Wavelength (nm).
wavelength = 488;
% Refractive index.
refInd = 1.33;

% === Reconstruct ===
% List of phase shifts.
phase = [0, -(2/5)*pi, -(4/5)*pi, -(6/5)*pi, -(8/5)*pi];
% Initial Kp values.
kpInit = [337, 400, 337, 463];
% Apodization radius (1/nm).
rApo = 336;
% Radius of the mask.
maskOD = 85.8;
maskID = 70.2;
% Pixel size of the PSF (nm).
pxSize = 102;

%% List the dataset.
if ~exist(inDir, 'dir')
    error(MSG_ID, 'Input directory does not exist.');
else
    l = dir(inDir);
    l = rmhiddendir(l);
    
    nData = numel(l);
    
    % Check if it's empty.
    if nData == 0
        warning(MSG_ID, 'No data in the input directory.');
        return;
    end
end

%% Prepare the output folder.
% Create the folder if it doesn't exist.
if ~exist(outDir, 'dir')
    mkdir(outDir);
else
    l = dir(outDir);
    l = rmhiddendir(l);
    
    % Check if it's empty.
    if numel(l) ~= 0
        warning(MSG_ID, 'Output directory is not empty.');
    end
end

% Create log file.
fid = fopen(fullfile(outDir, 'recon.log'), 'w');

%% Register the cleanup function.
% Must close the file on exit.
cleanObj = onCleanup(@()fclose(fid));

%% Load the PSF. 
psf = single(imread(psfFilename));
psf = centerpsf(psf);

%% Retrieve the geometry of the file.
% Calibrated Kp value.
kpCal = [];
% Parse the width and height info using the first raw data.
info = imfinfo(filepath(inDir, posInit, 0, 1));
% Raw image dimension, [width, height].
rawDim = [info.Width, info.Height];

%% Generate apodization function.
Iapo = cosapo(rawDim, rApo);

%% Create the mask for cross correlation.
Imask = crossmask(rawDim*2-1);

%% Generate phase coefficient matrix.
% NOTE: Aphase can be improved with the precise intensity of each beam.
% A x = b
Aphase = zeros(illMax);
for i = 1:illMax
    Aphase(:, i) = [3, ...
                    exp(-1i * phase(i) * 2*pi), ...
                    exp( 1i * phase(i) * 2*pi), ...
                    exp(-1i * phase(i) *   pi), ...
                    exp( 1i * phase(i) *   pi)];
end
Aphase = Aphase / 9;

%% Start processing through files.
for fIdx = 1:nData
    zPos = posInit + posDelta*fIdx;
    % Iterate through the orientations, index starts from 0 instead of 1.
    for oriIdx = 0:oriMax-1
        % Summed image, wiped out the illumination pattern.
        Isum = zeros(rawDim);
        % Raw images, illumination differences are stored plane by plane.
        Iraw = zeros([rawDim, illMax], 'single');
        
        % Load the patterns.
        for illIdx = 1:illMax
            fPath = filepath(inDir, zPos, oriIdx, illIdx);
            Iraw(:, :, illIdx) = single(imread(fPath));
        end
        % NOTE: Why divide by 4?
        Iraw = Iraw / 4;
        % Sum up the raw data to provide wide field result.
        Isum = sum(Iraw, 3);
        
        % 1st deconvolution.
        for illIdx = 1:illMax
            Iraw(:, :, illIdx) = deconvlucy(Iraw(:, :, illIdx), psf, 10);
        end
        
        % FT and retrieve the domains.
        
    end
end

end

function list = rmhiddendir(list)
%RMHIDDENDIR Remove directories that start with '.'

i = 1;
while i <= numel(list)
    n = list(i).name;
    if n(1) == '.'
        % Remove current item and stay here.
        list(i) = [];
    else
        % Increment to next item.
        i = i+1;
    end
end

end

function psf = centerpsf(psf)
%CENTERPSF Center the PSF.

% Find the center.
[~, idx] = max(psf(:));
[x, y] = ind2sub(size(psf), idx);
% Crop the image to the center.
psf = psf(1:2*(x-1), 1:2*(y-1));
% Remove the background (1 S.D.).
me = mean(psf(:));
sd = std(psf(:));
psf = abs(psf - (me-sd));
% Scale the PSF to have an integral value of one.
psf = psf / sum(psf(:));

end

function pstr = filepath(dstr, pos, ori, ill)
%FILEPATH Retrieve the filepath in the data set according to the
%measurement position and target orientation.

posStr = sprintf('%.3f', pos);
illPat = ['%', num2str(ori), '4d.tif'];
illStr = sprintf(illPat, ill);

pstr = fullfile(dstr, posStr, illStr);

end

function [vx, vy] = origingrid(sz)
%ORIGINGRID Generate a meshgrid, centered at zero.

szHalf = floor(sz/2);
% Shift the grid range to center at zero instead of sz/2.
%   i) Odd: -floor(sz/2):floor(sz/2)
%   ii) Even: -floor(sz/2):floor(sz/2)-1
szMin = -szHalf;
szMax = szHalf - xor(mod(sz, 2), 1);
% Meshgrid correspond to the provided dimension.
[vx, vy] = meshgrid(szMin(1):szMax(1), szMin(2):szMax(2));

end

function J = cosapo(sz, r)
%COSAPO Generates a cosine apodization function.

[vx, vy] = origingrid(sz);

J = cos((pi/2) * (sqrt( vx.^2 + vy.^2 )/r));
% Clamp the value to [0, +inf).
J(J < 0) = 0;

end

function M = crossmask(sz, od, id)
%CROSSMASK Generate the annular mask for cross-correlation.

% Calculate squared radius from inner/outer diameter.
or = (od/2)^2;
ir = (id/2)^2;

[vx, vy] = origingrid(sz);
% Calculate the squared distance for each cell.
R2 = vx.^2 + vy.^2;
M = (R2 > ir) & (R2 < or);

end
