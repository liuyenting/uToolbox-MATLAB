function batchshear_sim
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

% SIM parameters
nOri = 1;
nPhase = 5;

%% ask for directories and image type
inDir = uigetdir('C:\', 'Where are the original data?');
% check whether the input directory exists
if exist(inDir, 'dir') ~= 7
    error('postproc:batchshear', ...
          'Input directory does not exist.');
end

% hetero images
answer = questdlg('Are the images all of same size?', ...
                  'Heterogeneous Data', ...
                  'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        isHetero = false;
    case 'No'
        isHetero = true;
end

outDir = uigetdir(inDir, 'Where to store the result?');
% check whether the output directory is empty
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
    % create the directory
    status = mkdir(outDir);
    if ~status
        error('postproc:batchsehar', ...
              'Unable to create the output directory.');
    end
end

%% load the file list
%TODO use datastore instead of barebone directory list
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
    
    oldImg = io.fread(inFilePath);
    
    %% split orientations and phases
    % split the orientations and phases
    oldImg = sim.opmajor(oldImg, nOri, nPhase);
    % temporary storage
    newImg = [];
    volSz = [];
    % iterate through each stacks
    for iOri = 1:nOri
        for iPhase = 1:nPhase
            %% processing
            if isHetero
                % renew the shearing object
                shObj = postproc.Shear();
                shObj.setacqparam(acqParam);
            end

            shObj.loadstack(oldImg(:, :, :, iPhase, iOri));
            shObj.execute();
            
            % allocate space if non-exist
            if isempty(newImg)
                T = shObj.retrieveresult();
                volSz = size(T);
                newImg = zeros([volSz, nPhase, nOri]);
                newImg(:, :, :, iPhase, iOri) = T; 
            else
                newImg(:, :, :, iPhase, iOri) = shObj.retrieveresult(); %#ok<AGROW>
            end
        end
    end
  
    %% save the result
    % desired output is uint16 
    newImg = uint16(newImg);
    
    % revert to 3-D
    newImg = permute(newImg, [1, 2, 4, 5, 3]);
    newImg = reshape(newImg, volSz(1), volSz(2), []);
    
    outFilePath = fullfile(outDir, fileName);
    io.fwrite(newImg, outFilePath, 'Overwrite', true);
    
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
