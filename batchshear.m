function batchshear
%BATCHSHEAR Perform shearing on all the stacks in assigned directory.
%   BATCHSHEAR will first ask for input and output directory and then
%   perform the operation.

prompt = {'Objective angle (degree)', ...
          'Z step size (um)', ...
          'Pixel size (um)'};
title = 'Acquisition Parameters';
defAns = {'32.8', '0.5', '0.102'};
answer = inputdlg(prompt, title, 1, defAns);

if isempty(answer)
    return;
end

% acquisition paramters
acqParam.ObjectiveAngle = str2double(answer{1});   % [deg]
acqParam.ZStepWidth = str2double(answer{2});       % [um]
acqParam.PixelWidth = str2double(answer{3});       % [um]

% hetero images
% Construct a questdlg with three options
answer = questdlg('Are the images all of same size?', ...
                  'Heterogeneous Data', ...
                  'Yes', 'No', 'Yes');
% Handle response
switch answer
    case 'Yes'
        isHetero = false;
    case 'No'
        isHetero = true;
end

%% verify inputs
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

%% load the file list
list = dir(fullfile(inDir, '*.tif*'));
list = {list.name};

%% iterate through files
% Ignore warnings for unknown tags.
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

if ~isHetero
    % create the shearing object for sharing
    shObj = postproc.Shear();
    shObj.setacqparam(acqParam);
end

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
    
    %% processing
    if isHetero
        % renew the shearing object everytime
        shObj = postproc.Shear();
        shObj.setacqparam(acqParam);
    end

    shObj.loadstack(oldImg);
    shObj.execute();
    newImg = shObj.retrieveresult();
    
    %% save the result
    % desired output is uint16 
    %TODO probe the input image type and cast to it
    newImg = uint16(newImg);
    
    outFilePath = fullfile(outDir, fileName);
    tiff.imsave(newImg, outFilePath, true);
    
    %TODO debug the reason why files are not properly close by +tiff (dep)
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
        
    %% free the resources
    if isHetero
        % force invoke the destructor
        delete(shObj)
    end
end
tElapse = toc(tOuter);
fprintf('%f seconds to process the incoming data\n', tElapse);

warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

end