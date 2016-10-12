function pprecon3
%PPRECON3 Reconstruct the 3D SIM dataset using plane-by-plane method.

clearvars; close all;

MSG_ID = 'simrecon:pprecon3';

%% Preset the parameters.
inDir = '~/Documents/SIM/RAWcell1';
psfFilename = '~/Documents/SIM/ygbead_zp1um_NAp55nap44_ExpPsf.tif';

[path, name, ~] = fileparts(inDir);
dateStr = datestr(now, 'yyyymmddHHMMss');
outDir = fullfile(path, [name, '_', dateStr]);

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
        warning(MSG_ID, 'Output directory is not empty.');
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
info = imfinfo(filepath(inDir, posInit, 0, 1));
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
    p = phase(i);
    Aphase(:, i) = [3, ...
                    exp(-1i * p * 2*pi), exp(1i * p * 2*pi), ...
                    exp(-1i * p * pi),  exp(1i * p * pi)];
end
Aphase = Aphase / 9;

%% Start processing through files.
% DEBUG
nData = 1;
for fIdx = 1:nData
    tic;
    
    zPos = posInit + (posDelta*1e-3)*fIdx;
    % Iterate through the orientations, index starts from 0 instead of 1.
    for oriIdx = 0:oriMax-1
        % Summed image, wiped out the illumination pattern.
        Isum = zeros(rawDim, 'single');
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
        Isum = mean(Iraw, 3);

        % 1st deconvolution.
        for illIdx = 1:illMax
            Iraw(:, :, illIdx) = deconvlucy(Iraw(:, :, illIdx), psf, 10);
        end

        % FT and retrieve the domains.
        E = zeros(size(Iraw), 'single');
        for illIdx = 1:illMax
            E(:, :, illIdx) = fftshift(fft2(Iraw(:, :, illIdx)));
            % Limit by the apodization function.
            E(:, :, illIdx) = E(:, :, illIdx) .* Iapo;
        end

        % Isolate the results.
        %TODO: Directly reshape from 3D to 2D.
        Eflat = zeros([nElem, illMax], 'single');
        for illIdx = 1:illMax
            Eflat(:, illIdx) = reshape(E(:, :, illIdx), [nElem, 1]);
        end
        Dflat = Aphase \ Eflat.';
        Dflat = Dflat.';

        % Resize back the frequency domains, and re-use E.
        %TODO: Directly reshape from 3D to 2D.
        D = zeros([rawDim, illMax], 'single');
        for illIdx = 1:illMax
            D(:, :, illIdx) = reshape(Dflat(:, illIdx), rawDim);
        end
        
        % Pad the images, only pad along the image planes.
        D = padarray(D, [rawDim/2, 0]);
    
        % Using cross-correlation to find the Kp through the m1 terms.
        %Icross = zeros([2*rawDim, ilMax], 'single');
        % TODO: Bypass.
        
        shift = zeros([2, 4], 'single');
        phaseInit = zeros([1, 4], 'single');
        
        % Restore the wide field image for comparison.
        Iwf = abs(ifft2(ifftshift(D(:, :, 1))));
        
        % Preset the figure for temporary result.
        f = figure('Name', '', 'NumberTitle', 'off');
        f.Position = [100, 100, 1024, 512];
%         f.MenuBar = 'none';
%         f.ToolBar = 'none';
        % Delta phase.
        dp = 2*pi / phStepRes;
        phaseOpt = zeros(2, phStepRes);
        for pIdx = 1:phStepRes
            s = sprintf('Localize @ %.3frad', pIdx*dp);
            f.Name = s;
            
            % Iterate through m2~m4 terms.
            sgn = 1;
            for m = 2:illMax
                % "Increment" the phase.
                %TODO: Apply at the shifted region instead of overall
                %region.
                D(:, :, m) = D(:, :, m) * exp(sgn * 1i * dp);
                % Toggle the sign for next cycle.
                sgn = -sgn;
            end
            
            % Test for m2/m3 terms.
            Itmp = D(:, :, 1) + D(:, :, 2) + D(:, :, 3);
            % Revert to spatial domain.
            Itmp = abs(ifft2(ifftshift(Itmp)));
            Itmp = Itmp .* Iwf;
            phaseOpt(1, pIdx) = sum(Itmp(:));
            
            % Test for m4/m5 terms;
            Itmp = D(:, :, 1) + D(:, :, 4) + D(:, :, 5);
            % Revert to spatial domain.   
            Itmp = abs(ifft2(ifftshift(Itmp)));
            Itmp = Itmp .* Iwf;
            phaseOpt(2, pIdx) = sum(Itmp(:));
        end
        % Move the data to the optimal phase.
        [~, maxIdx] = max(phaseOpt, [], 2);
        dp = (phStepRes-maxIdx) * dp;
        % Apply through m2~m4 terms.
        sgn = -1;
        for m = 2:illMax
            % "Increment" the phase.
            D(:, :, m) = D(:, :, m) * exp(sgn * 1i * dp);
            % Toggle the sign for next cycle.
            sgn = -sgn;
        end
        
        % Create SIM reconstructed PSF.
        
        % Create the pattern and their respective raw data.
        
        % Reconstruct the result.
        
        % Find the initial phases for m_i terms.
        
        % Save the frequency domain results.
        
        % Reconstruct the PSF.
        
        % Export the final image.
        
        % Export SIM image with all the orientations.
        
        % Second deconvolution.
        
        % Write back the Kp.
        
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
