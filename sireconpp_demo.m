%SIRECONPP_DEMO Demonstrate the plan-by-plan reconstruction feature.
%
%   TBA

% close all the figures
close all;
% clear variables from the workspace
clearvars; 
% wipe persistent variables
clearvars -global;

%TODO remove OS dependent test code
if ispc
    srcDir = 'D:\\Andy\\07252017_SIMforAndy\\cell3_3DSIM_obj';
    psfFileName = 'D:\\Andy\\07252017_SIMforAndy\\PSF3_3DSIM_obj\\RAWpsf3a3DSIM_ch0_561nm_cropped.tif';
else
    srcDir = 'data/07252017_SIMforAndy/cell3_3DSIM_obj';
    psfFileName = 'data/07252017_SIMforAndy/PSF3_3DSIM_obj/RAWpsf3a3DSIM_ch0_561nm_cropped.tif';

    % prepend the root
    srcDir = fullfile(userpath, srcDir);
    psfFileName = fullfile(userpath, psfFileName);
end

dstSuffix = '_sirecon';
dstDir = [srcDir, dstSuffix];

%% setup parameters
siparms = struct;

siparms.Debug = true;
siparms.DebugPath = fullfile(dstDir, 'debug');

% acquisition info
siparms.Orientations = 1;
siparms.Phases = 5;

% beam intensity
siparms.I0 = 1;
siparms.I1 = 1;

% transfer function
siparms.PixelSize = 100;            % [nm] or [nm, nm]

siparms.RefractiveIndex = 1.33;     % refractive index of the medium
siparms.Wavelength = 520;           % [nm]
siparms.NA = 1.1;                   % effective numerical aperture

% deconvolution
siparms.PreDeconv = 5;
siparms.PostDeconv = 5;
siparms.PadSize = 10;

% phase retrieval
siparms.ApodizeRatio = 0.5;
siparms.RetrievalInterpRatio = 2;

% generalized Wiener filter
siparms.WienerConstant = 1e-3;

%% verify the input
% check whether the input directory exsists
if exist(srcDir, 'dir') ~= 7
    error(generatemsgid('InvalidInDir'), ...
          'Input directory does not exist.');
end

% check whether the output directory is empty
if exist(dstDir, 'dir') == 7
    content = dir(dstDir);
    if numel(content) > 2
        if siparms.Debug
            warning(generatemsgid('InvalidOutDir'), ...
                    'Output directory is not empty.');
            % recreate the folder
            util.rmcontent(dstDir);
        else
            error(generatemsgid('InvalidOutDir'), ...
                  'Output directory is not empty.');
        end
    end
end

% duplicate the pixel size if isotropic
if length(siparms.PixelSize) == 1
    siparms.PixelSize = siparms.PixelSize * [1, 1];
end

%% create the destination directory
% destination root
status = mkdir(dstDir);
if ~status
    error(generatemsgid('InvalidOutDir'), ...
          'Unable to create the output directory.');
end

% create the directory
if siparms.Debug
    mkdir(siparms.DebugPath);
end

%% load the file list
%TODO load 488nm only for now
imds = imageDatastore(fullfile(srcDir, '*c3a*.tif*'), ...
                      'ReadFcn', @tiff.TIFFStack);

%% preload the PSF
nOri = siparms.Orientations;
nPhases = siparms.Phases;

Ipsf = tiff.imread(psfFileName);
Ipsf = single(Ipsf);

% re-order the stack to orientation-wise and phase-wise
Ipsf = sim.opmajor(Ipsf, nOri, nPhases);

psfSz = size(Ipsf);
psfnz = psfSz(3);
% select the mid-plane for 5P-2D-SIM
Ipsf = Ipsf(:, :, floor(psfnz/2)+1, :, :);
% remove the Z dimension
psfSz = psfSz(1:2);
Ipsf = reshape(Ipsf, [psfSz, nPhases, nOri]);

% center the PSF
for iOri = 1:nOri
    for iPhase = 1:nPhases
        Ipsf(:, :, iPhase, iOri) = image.centerpsf(Ipsf(:, :, iPhase, iOri));
    end
end

if siparms.Debug
    figure('Name', 'PSF', 'NumberTitle', 'off');
    for iOri = 1:nOri
        for iPhase = 1:nPhases
            subplot(nOri, nPhases, (iOri-1)*nOri+iPhase);
            imagesc(Ipsf(:, :, iPhase, iOri));
                axis image;
                title(sprintf('O%d P%d', iOri, iPhase));
        end
    end
end

% save into SI parameter
siparms.PSF = Ipsf;
siparms.TransFunc = [];

%TODO identify 2-D or 3-D SIM from the dimension of PSF

%% iterate through the files
% ignore warnings for unknown tags for current session
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

profile on;
profile off;
   
tOuter = tic;
while hasdata(imds)
    tInner = tic;
    fprintf('\n');
    
    % load the file
    [I, info] = read(imds);
    fPath = info.Filename;
    
    [~, fName, fExt] = fileparts(fPath); 
    fprintf('[%s]\n', fName);
    
    % convert to floating point
    I = single(I);
    % re-order the stack to orientation-wise and phase-wise
    I = sim.opmajor(I, siparms.Orientations, siparms.Phases);
    
    % normalize the intensity across phases
    I = sim.normopint(I, siparms);
    % execute
    J = sim.sireconpp(I, siparms);
    
    % save the volume
    fPath = fullfile(dstDir, [fName, '.', fExt]);
    tiff.imsave(J, fPath);
    
    % cleanup
    fclose('all');
    
    tElapse = toc(tInner);
    fprintf('%.3fs elapsed\n', tElapse);
end
tElapse = toc(tOuter);
fprintf('+++ %.3fs to complete the time lapse +++\n', tElapse);

% profile viewer;