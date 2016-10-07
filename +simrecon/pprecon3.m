function pprecon3
%PPRECON3 Reconstruct the 3D SIM dataset using plane-by-plane method.

clearvars; close all;

MSG_ID = 'simrecon:pprecon3';

%% Preset the parameters.
inDir = '~/Documents/SIM/RAWcell1';
psfFilename = '~/Documents/SIM/ygbead_zp1um_NAp55nap44_ExpPsf.tif';

[path, name, ~] = fileparts(inDir);
outDir = fullfile(path, [name, '_recon']);

% Wavelength (nm).
wavelength = 488;
% Refractive index.
refInd = 1.33;

% List of phase shifts.
phase = [0, -(2/5)*pi, -(4/5)*pi, -(6/5)*pi, -(8/5)*pi];
% Initial Kp values.
kpInit = [337, 400, 337, 463];
% Apodization radius (1/nm).
rApo = 336;
% Radius of the mask.
rMask = 39;
% Pixel size of the PSF (nm).
pxSize = 102;

% Maximum of orientation.
oriMax = 5;

%% List the dataset.
if ~exist(inDir, 'dir')
    error(MSG_ID, 'Input directory does not exist.');
else
    dirList = dir(inDir);
    dirList = rmhiddendir(dirList);
    
    % Check if it's empty.
    if numel(dirList) == 0
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
cleanObj = onCleanup(@()fclose(fid));

%% Load the PSF. 
psf = single(imread(psfFilename));
psf = centerpsf(psf);

end

%% Retrieve the geometry of the file.
% Calibrated Kp value.
kpCal = [];
% Dummy read the file.

eField = single(tiff.imread());

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