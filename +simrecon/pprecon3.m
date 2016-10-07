function pprecon3
%PPRECON3 Reconstruct the 3D SIM dataset using plane-by-plane method.

clearvars; close all;

MSG_ID = 'simrecon:pprecon3';

%% Preset the parameters.
inDir = 'RAWcell1';
psfFilename = 'ygbead_zp1um_NAp55nap44_ExpPsf.tif';

[path, name, ~] = fileparts(inDir);
outDir = fullfile(path, [name, '_recon']);

% Wavelength (nm).
wavelength = 488;
% Refractive index.
refInd = 1.33;

% List of phase shifts.
phase = [0, -(2/5)*pi, -(4/5)*pi, -(6/5)*pi, -(8/5)*pi];
% Kp values.
kp = [337, 400, 337, 463];
% Apodization radius (1/nm).
rApo = 336;
% Radius of the mask.
rMask = 39;
% Pixel size of the PSF (nm).
pxSize = 102;

%% Prepare the output folder.
% Create the folder if it doesn't exist.
if ~exist(outDir, 'dir')
    mkdir(outDir);
else
    dirList = dir(outDir);
    dirList = rmhiddendir(dirList);
    
    % Check if it's empty.
    if numel(dirList) ~= 0
        error(MSG_ID, 'Output directory is not empty.');
    end
end

% Create log file.


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
