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
    srcDir = 'D:\\Andy\\06302017_SIM\\Tan1_3_3DSIM';
    psfFileName = 'D:\\Andy\\06302017_SIM\\PSF\\PSF_NA1p1_520nm_x100nm.tif';
else
    srcDir = 'data/sim/06302017_SIM/Tan1_3_3DSIM';
    psfFileName = 'data/sim/06302017_SIM/PSF/PSF_NA1p1_520nm_x100nm.tif';

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

% kp FFT upsampling ratio
siparms.KpUpsamplingRatio = 2;

% deconvolution
siparms.PreDeconv = 5;
siparms.PostDeconv = 5;
siparms.PadSize = 10;

% phase retrieval
siparms.ApodizeRatio = 0.5;
siparms.RetrievalInterpRatio = 2;

% generalized Wiener filter
siparms.NA = 1.1;           % effective NA
siparms.Wavelength = 520;   % [nm]
siparms.PixelSize = 100;    % [um] or [um, um]
siparms.WienerConstant = 0.001;

%% verify the input
% check whether the input directory exsists
if exist(srcDir, 'dir') ~= 7
    error(generatemsgid('InvalidInput'), ...
          'Input directory does not exist.');
end

% check whether the output directory is empty
if exist(dstDir, 'dir') == 7
    content = dir(dstDir);
    if numel(content) > 2
        if siparms.Debug
            warning(generatemsgid('InvalidOutput'), ...
                    'Output directory is not empty.');
            % recreate the folder
            util.rmcontent(dstDir);
        else
            error(generatemsgid('InvalidOutput'), ...
                  'Output directory is not empty.');
        end
    end
end

% duplicate pixel size if isotropic
if length(siparms.PixelSize) == 1
    siparms.PixelSize = siparms.PixelSize * [1, 1];
end

%% create the destination directory
% destination root
status = mkdir(dstDir);
if ~status
    error(generatemsgid('InvalidOutput'), ...
          'Unable to create the output directory.');
end

% create the directory
if siparms.Debug
    mkdir(siparms.DebugPath);
end

%% load the file list
%TODO load 488nm only for now
imds = imageDatastore( ...
    fullfile(srcDir, '*488*.tif*'), ...
    'ReadFcn', @tiff.TIFFStack ...
);

%% preload the PSF
Ipsf = tiff.imread(psfFileName);
Ipsf = single(Ipsf);

% center the PSF
Ipsf = image.centerpsf(Ipsf);

% normalize PSF, ensure kernel sum is 1
Ipsf = Ipsf / sum(Ipsf(:));

% save into SI parameter
siparms.PSF = Ipsf;

%% iterate through the files
% ignore warnings for unknown tags for current session
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

profile on;
profile off;
   
tOuter = tic;
while hasdata(imds)
    tInner = tic;
    
    % load the file
    [I, info] = read(imds);
    imSz = size(I);
    fPath = info.Filename;
    
    [~, fName, fExt] = fileparts(fPath); 
    fprintf('[%s]\n', fName);
    
    % convert to floating point
    I = single(I);
    % re-order the stack to phase-wise
    I = sim.opmajor(I, siparms.Orientations, siparms.Phases);
    
    % normalize the intensity across phases
    I = sim.normopint(I, siparms);
    % execute
    J = sim.sireconpp(I, volSz, siparms);
    
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