clearvars; close all;

%% user configurations

% directory that stored the raw images, absolute path
inputDir = uigetdir('Where is the source directory?');

% suffix for the output directory
%   If the postfix is "_deskew" and "C:\Data" as the input directory,
%   "C:\Data_deskew" will be the output directory.
outputDirPostfix = '_deskew';

% acquisition parameters
pixelWidth = 0.102; % [um]
zStepWidth = 0.5; % [um]
objectiveAngle = 32.8; % [deg]

%% precalculate the translation matrix and result image size
objectiveAngle = degtorad(objectiveAngle);
%objectiveAngle = objectiveAngle * pi/180;

% d_xy = d_z * sin( pi/2 - theta )
xyUnitShift = zStepWidth * cos(objectiveAngle) / pixelWidth;
clear pixelWidth zStepWidth objectiveAngle;

%% Acquire file list and create output directory.
fileList = dir(fullfile(inputDir, '*.tif*'));
totalFiles = numel(fileList);

% generate output file path
outputDir = [inputDir, outputDirPostfix];
[~, ~, messageId] = mkdir(outputDir);
if strcmp(messageId, 'MATLAB:MKDIR:DirectoryExists')
    warning(generatemsgid('OutputDirectoryExists'), ...
            'Please remove the directory "%s".', outputDir);
    return;
end
clear messageId outputDirPostfix;

%% Iterate through the stacks and process them.
isFirstRun = true;
outImgStack = [];

warning('off', 'all');
tic;
for fileIdx = 1:totalFiles
    file = fileList(fileIdx);
    
    % read the layers in selected stack
    fullInFilePath = fullfile(inputDir, file.name);
    tiffStack = io.fread(fullInFilePath);
    
    stackDim = size(tiffStack);
    totalLayers = stackDim(3);
    
    if isFirstRun
        isFirstRun = false;
        
        % generate reference frame dimension
        height = int32(stackDim(1));
        width = int32(stackDim(2) + xyUnitShift*totalLayers);
        layerDim = [height, width];
        imrefCls = imref2d(layerDim, ...
                           [0, single(width)], [0, single(height)]);
        
        % pre-allocate
        outImgStack = zeros([height, width, totalLayers]);
                       
        clear width height layerDim;
    end
    
    % create output file name
    fullOutFilePath = fullfile(outputDir, file.name);
    
    transMat = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    for idx = 1:totalLayers
        % create the translation description object
        transObj = affine2d(transMat);

        % translate the image
        outImg = imwarp(tiffStack(:,:,idx), transObj, 'linear', ...
                        'OutputView', imrefCls, 'FillValues', 0);
                    
        % write into the stack
        outImgStack(:, :, idx) = outImg; %#ok<SAGROW>

        % update translation matrix
        transMat = transMat + [0, 0, 0; 0, 0, 0; xyUnitShift, 0, 0];
    end
    
    % save as TIFF stack
    outImgStack = uint16(outImgStack);
    io.fwrite(outImgStack, fullOutFilePath);

    % update the progress
    fprintf('%d/%d\n', fileIdx, totalFiles);
end
toc;
warning('on', 'all');

clearvars; 