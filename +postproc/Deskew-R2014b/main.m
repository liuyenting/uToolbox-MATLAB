clearvars; close all;

%% User configurations

% Directory that stored the raw images, absolute path.
%inputDir = '/Users/Andy/Documents/Deskew-Resources/Data/CRY2_clone4-4_zp5um_20ms_a2s_r2s_Raw';
inputDir = 'F:\07122016_CHOK1_WITAT\cell4_zp5um_20ms_a9p5_rp5s';

% Postfix for the output directory.
%   If the postfix is "_deskew" and "C:\Data" as the input directory,
%   "C:\Data_deskew" will be the output directory.
outputDirPostfix = '_deskew';

% Background noise characteristics.
useNoise = true;
noiseAverage = 90.4633;
noiseSD = 0.3053;
% If useNoise is false, background will be filled by this scalar.
fillScalar = 42;

% Acquisition parameters.
pixelWidth = 0.102; % [um]
zStepWidth = 0.5; % [um]
objectiveAngle = 32.8; % [deg]

% Enable CUDA.
useCuda = false;

%% Add the libraries.
% Library paths will only persist in current MATLAB session.
libPath = genpath('lib');
addpath(libPath);
clear libPath;

%% Precalculate the translation matrix and result image size.
objectiveAngle = degtorad(objectiveAngle);
%objectiveAngle = objectiveAngle * pi/180;

% d_xy = d_z * sin( pi/2 - theta )
xyUnitShift = zStepWidth * cos(objectiveAngle) / pixelWidth;
clear pixelWidth zStepWidth objectiveAngle;

%% Acquire file list and create output directory.
fileList = dir(fullfile(inputDir, '*.tif*'));
totalFiles = numel(fileList);

% Generate output file path.
outputDir = [inputDir, outputDirPostfix];
[~, ~, messageId] = mkdir(outputDir);
if strcmp(messageId, 'MATLAB:MKDIR:DirectoryExists')
    warning('Deskew:OutputDirectoryExists', ...
            'Please remove the directory "%s".', outputDir);
    return;
end
clear messageId outputDirPostfix;

%% Iterate through the stacks and process them.
isFirstRun = true;
noiseTemplate = []; 

% TIFF file options.
tiffOption.append = true;

% Create progress bar.
progressbar('Files', 'Layers');

warning('off', 'all');
tic;
for fileIdx = 1:totalFiles
    file = fileList(fileIdx);
    
    % Read the layers in selected stack.
    fullInFilePath = fullfile(inputDir, file.name);
    tiffStack = TIFFStack(fullInFilePath);
    
    stackDim = size(tiffStack);
    totalLayers = stackDim(3);
    
    if isFirstRun
        isFirstRun = false;
        
        % Generate reference frame dimension.
        height = int32(stackDim(1));
        width = int32(stackDim(2) + xyUnitShift*totalLayers);
        layerDim = [height, width];
        imrefCls = imref2d(layerDim, ...
                           [0, single(width)], [0, single(height)]);

        if useNoise
            disp('Generating noise...');
            
            % Calculate noise SNR.
            snr = noiseAverage/noiseSD;
            snr_db = mag2db(snr);
            
            % Generate noise template.
            noiseTemplate = noiseAverage*ones(1, width*height);
            noiseTemplate = simulate.awgnoise(noiseTemplate, snr_db);
            
            % Reshape to image size.
            noiseTemplate = reshape(noiseTemplate, layerDim);
            
            % Set template fill value.
            fillScalar = NaN;
            
            clear snr snr_db;
            clear noiseAverage noiseSD;
        end
        
        clear width height layerDim;
    end
    
    % Create output file name.
    fullOutFilePath = fullfile(outputDir, file.name);
    
    transMat = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    for idx = 1:totalLayers
        % Create the translation description object.
        transObj = affine2d(transMat);

        % Translate the image.
        outImg = imwarp(tiffStack(:,:,idx), transObj, 'linear', ...
                        'OutputView', imrefCls, 'FillValues', fillScalar);
        if useNoise
            % Fill in the noise.
            fillLocation = ~outImg;
            outImg(fillLocation) = noiseTemplate(fillLocation);
        end
        
        % Save to the TIFF stack.
        saveastiff(outImg, fullOutFilePath, tiffOption);
        
        % Update layer progress.
        progressbar([], idx/totalLayers);
        
        % Update translation matrix.
        transMat = transMat + [0, 0, 0; 0, 0, 0; xyUnitShift, 0, 0];
    end
    
    % Update file progress.
    progressbar(fileIdx/totalFiles);
end
toc;
warning('on', 'all');
