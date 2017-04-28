function batchshear
%BATCHSHEAR Perform shearing on all the stacks in assigned directory.
%   BATCHSHEAR will first ask for input and output directory and then
%   perform the operation.

% Acquisition paramters.
acqParam.ObjectiveAngle = -32.8; % [deg]
acqParam.ZStepWidth = 0.3;      % [um]
acqParam.PixelWidth = 0.102;    % [um]

%% Verify input.
% Ask for directory.
inDir = uigetdir('C:\', 'Where are the original data?');
outDir = uigetdir(inDir, 'Where to store the result?');

% Check whether the input directory exists.
if exist(inDir, 'dir') ~= 7
    error('postproc:batchshear', ...
          'Input directory does not exist.');
end

% Check whether the output directory is empty.
if exist(outDir, 'dir') == 7
    content = dir(outDir);
    if numel(content) > 2
        button = questdlg('Output directory is not empty.', ...
                          'Continue the procedure?', ...
                          'Continue', 'Abort', 'Abort');
        if strcmp(button, 'Abort')
            error('postproc:batchshear', ...
                  'Output directory is not empty.');
        end
    end
else
    % Create the directory.
    status = mkdir(outDir);
    if ~status
        error('postproc:batchsehar', ...
              'Unable to create the output directory.');
    end
end

%% Load the file list.
list = dir(fullfile(inDir, '*.tif*'));
list = {list.name};

%% Iterate through files.
% Ignore warnings for unknown tags.
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

nFile = numel(list);
tOuter = tic; tAverage = 0;
for iFile = 1:nFile
    tInner = tic;
    
    fileName = list{iFile};
    inFilePath = fullfile(inDir, fileName);
    
    oldImgStack = tiff.TIFFStack(inFilePath);
    [nx, ny, nz] = size(oldImgStack);
    oldImg = oldImgStack(:);
    oldImg = reshape(oldImg, nx, ny, nz);
    
    %% Begin processing.
    % Create the shearing object.
    shObj = postproc.Shear();
    shObj.setacqparam(acqParam);
    % Load and run.
    shObj.loadstack(oldImg);
    shObj.execute();
    newImg = shObj.retrieveresult();
    
    %% Save the result.
    newImg = uint16(newImg);
    
    outFilePath = fullfile(outDir, fileName);
    tiff.imsave(newImg, outFilePath, true);
    
    fclose('all');
    
    tElapse = toc(tInner);
    if tAverage == 0
        tAverage = tElapse;
    else
        tAverage = (tAverage+tElapse)/2;
    end
    tEta = (nFile-iFile) * tAverage;
    strEta = datestr(tEta/86400, 'HH:MM:SS');
    fprintf('%d file processed, %f seconds elapsed, ETA %s\n', ...
            iFile, tElapse, strEta);
end
tElapse = toc(tOuter);
fprintf('%f seconds to process the incoming data\n', tElapse);

warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

end