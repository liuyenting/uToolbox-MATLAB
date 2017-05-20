function pprecon3
%PPRECON3 Reconstruct the 3D SIM dataset using plane-by-plane method.

clearvars; close all;

global MSG_ID;
MSG_ID = 'simrecon:pprecon3';

%% Preset the parameters.
%inDir = '~/Documents/SIM/RAWcell1';
%psfFilename = '~/Documents/SIM/ygbead_zp1um_NAp55nap44_ExpPsf.tif';
inDir = 'C:\Users\Avizo\Downloads\sim_test_data\RAWcell1';
psfFilename = 'C:\Users\Avizo\Downloads\sim_test_data\ygbead_zp1um_NAp55nap44_ExpPsf.tif';

[path, name, ~] = fileparts(inDir);
dateStr = datestr(now, 'yyyymmddHHMMss');
outDir = fullfile(path, [name, '_', dateStr]);

% === Data Sets ===
% Initial position (um).
posInit = 100;
% Position increment (nm).
posDelta = 250;
% Maximum of orientation types.
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
phaseCoeff = [0, -2/5, -4/5, -6/5, -8/5];
% Initial Kp values, each row represents an orientation.
kpInit = [337, 463, 337, 400]; % Inverted.
% Apodization radius (1/nm).
rApo = 336;
% Radius of the mask.
maskOD = 85.8;
maskID = 70.2;
% Resolution of phase steps (#).
phStepRes = 20;
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
        error(MSG_ID, 'Output directory is not empty.');
    end
end

% Create log file.
fid = fopen(fullfile(outDir, 'recon.log'), 'w');

%% Register the cleanup function.
% Must close the file on exit.
% TODO: Modify to rename the directory from UUID to datetime-based.
cleanObj = onCleanup(@()fclose(fid));

%% Load the PSF.
psf = single(imread(psfFilename));
psf = centerpsf(psf);

% f = figure('Name', 'PSF', 'NumberTitle', 'off');
% f.Position = [100, 100, 3*size(psf)];
% f.MenuBar = 'none';
% f.ToolBar = 'none';
% subplot('Position', [0 0 1 1]);
%     imagesc(psf);
%     axis equal off;

%% Retrieve the geometry of the file.
% Calibrated Kp value.
kpCal = [];
% Parse the width and height info using the first raw data.
info = imfinfo(filepath(inDir, posInit, 1, 1));
% Raw image dimension, [width, height].
rawDim = [info.Width, info.Height];
% N elements in a single plane.
nElem = prod(rawDim);

%% Generate apodization function.
Iapo = cosapo(rawDim, rApo);

%% Create the mask for cross correlation.
Imask = crossmask(rawDim*2-1, maskOD, maskID);

%% Generate phase coefficient matrix.
% NOTE: Aphase can be improved with the precise intensity of each beam.
% A x = b
Aphase = zeros(illMax);
for i = 1:illMax
    p = phaseCoeff(i);
    Aphase(i, :) = [3, ...
                    exp(-1i * p * 2*pi), exp(1i * p * 2*pi), ...
                    exp(-1i * p * pi),  exp(1i * p * pi)];
end
% TODO: Why divide by 9.
Aphase = Aphase / 9;

%% Start processing through files.
% DEBUG
nData = 1;
for fIdx = 1:nData
    tic;
    
    zPos = posInit + (posDelta*1e-3)*fIdx;
    % Iterate through the orientations, index starts from 0 instead of 1.
    for oriIdx = 1:oriMax
        % Summed image, wiped out the illumination pattern.
        Isum = zeros(rawDim);
        % Raw images, illumination differences are stored plane by plane.
        Iraw = zeros([rawDim, illMax]);

        %% Load the patterns.
        for illIdx = 1:illMax
            fPath = filepath(inDir, zPos, oriIdx, illIdx);
            Iraw(:, :, illIdx) = single(imread(fPath));
        end
        %TODO: Scale to [0, 1] and remember the min/max.
        Iraw = Iraw / 4;
        % Sum up the raw data to provide wide field result.
        %TODO: Move out the mean procedure.
        Isum = mean(Iraw, 3);

        %% 1st deconvolution.
        for illIdx = 1:illMax
            Iraw(:, :, illIdx) = deconvlucy(Iraw(:, :, illIdx), psf, 10);
        end

        %% FT and retrieve the domains.
        E = zeros(size(Iraw));
        for illIdx = 1:illMax
            E(:, :, illIdx) = fftshift(fft2(Iraw(:, :, illIdx)));
            % Limit by the apodization function.
            E(:, :, illIdx) = E(:, :, illIdx) .* Iapo;
        end

        %% Solve by the coefficients.
        E = reshape(E, [nElem, illMax]);
        D = Aphase \ E.';
        D = D.';
        D = reshape(D, size(Iraw));
        
        % Pad the images, only pad along the image planes.
        S = padarray(D, [rawDim/2, 0]);
        
        %% Restore the wide field image for comparison.
        Iwf = abs(ifft2(ifftshift(S(:, :, 1))));
        figure('Name', 'Wide-field', 'NumberTitle', 'off');
            imagesc(Iwf);
            axis equal tight;
        % Flatten the WF image, so we won't tamper with Itmp later on.
        Iwf = Iwf(:);
        
        %% Calculate the shift positions for the terms.
        % TODO: Calibrate the kp and generate kpCal.
        mShift = zeros(illMax, 2);
        for i = 2:2:illMax
            mShift(i, :) = [kpInit(oriIdx, i/2+1)-1, kpInit(oriIdx, i/2)-1];
            mShift(i+1, :) = rawDim - mShift(i, :);
        end
        kpCal = kpInit;

        %% Search for the inital phase offsets.
        % List of phases to test for.
        phaseOff = linspace(0, 2*pi, phStepRes);
        % Optimal phased image.
        Iopt = zeros(rawDim*2);
        % Plus/Minus temporary phased image.
        Sp = zeros(rawDim*2);
        Sm = zeros(rawDim*2);
        for m = 2:2:illMax
            maxCorrVal = 0;
            maxCorrPh = 0;
            for p = phaseOff
                Sp(mShift(m, 1)+1:mShift(m, 1)+rawDim(1), ...
                   mShift(m, 2)+1:mShift(m, 2)+rawDim(2)) ...
                    = D(:, :, m) * exp(1i * p);
                Sm(mShift(m+1, 1)+1:mShift(m+1, 1)+rawDim(1), ...
                   mShift(m+1, 2)+1:mShift(m+1, 2)+rawDim(2))...
                    = D(:, :, m+1) * exp(-1i * p);
                
                % Test for correlation result, update if overlapped more.
                Stmp = S(:, :, 1) + Sp + Sm;
                Itmp = abs(ifft2(ifftshift(Stmp)));
                s = sum(Itmp(:) .* Iwf);
                % Update the maximum result if greater.
                if s > maxCorrVal
                    % Update the maximum thresholds.
                    maxCorrVal = s;
                    maxCorrPh = p;
                    % Update the m-terms.
                    S(:, :, m) = Sp;
                    S(:, :, m+1) = Sm;
                    % Save the latest optimal result.
                    Iopt = Itmp;
                end
            end
            
            ts = sprintf('phase = %.3f\ncorr = %.5e', maxCorrPh, maxCorrVal);
            figure('Name', ['m', num2str(m/2), ' term'], 'NumberTitle', 'off');
                imagesc(Iopt);
                axis equal tight;
                title(ts, 'FontSize', 14);
        end

        %% Create the patterns.
        k = 2*pi / (wavelength/refInd);
        
        %% Reconstruct the result.
        
        %% Export SIM image with all the orientations.
        
        %% Create reconstructed PSF.
        
        %% Second deconvolution.
        
        %% Export the final image.
        
        %% Write back the Kp.
        
    end
    
    t = toc;
    s = sprintf('%.3fs to process layer %d', t, fIdx);
    disp(s);
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
% Remove the background (1 adjusted S.D.).
me = mean(psf(:));
sd = estnoise(psf);
psf = abs(psf - (me-sd));
% Scale the PSF to have an integral value of one.
psf = psf / sum(psf(:));

end

function s = estnoise(I)
%ESTNOISE Estimate the background noise level in terms of sigma.
%   J. Immerkær, Fast Noise Variance Estimation
%   Computer Vision and Image Understanding, Vol. 64, No. 2, Pg. 300-302

[h, w] = size(I);

% Compute sum of absolute values of Laplacian.
k = [1, -2, 1; -2, 4, -2; 1, -2, 1];
I = abs(conv2(I, k));
s = sum(I(:));

% Scale the sigma with proposed coefficients.
s = s * sqrt(pi/2) / (6*(w-2)*(h-2));

end

function pstr = filepath(dstr, pos, ori, ill)
%FILEPATH Retrieve the filepath in the data set according to the
%measurement position and target orientation.

global MSG_ID;

assert(ori > 0, MSG_ID, 'Orientation ID should be positive.');

posStr = sprintf('%.3f', pos);
illPat = ['%', num2str(ori-1), '4d.tif'];
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
