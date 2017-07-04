%SIRECONPP_DEMO Demonstrate the plan-by-plan reconstruction feature.
%
%   TBA

close all;
clearvars;

srcDir = 'F:\\Andy\\06302017_SIM\\Tan1_3_3DSIM';
psfFileName = 'F:\\Andy\\06302017_SIM\\PSF\\PSF_NA1p1_520nm_x100nm.tif';

dstSuffix = '_sirecon';
dstDir = [srcDir, dstSuffix];

% orientations
nOri = 1;
% phases
nPhase = 5;

%% verify the input
% check whether the input directory exsists
if exist(srcDir, 'dir') ~= 7
    error('sireconpp_demo', 'Input directory does not exist.');
end

% check whether the output directory is empty
if exist(dstDir, 'dir') == 7
    content = dir(dstDir);
    if numel(content) > 2
        error('sireconpp_demo', 'Output directory is not empty.');
    end
else
    % create the directory
    status = mkdir(dstDir);
    if ~status
        error('sireconpp_demo', 'Unable to create the output directory.');
    end
end

%% load the file list
%TODO only load 488nm for simplicity
list = dir(fullfile(srcDir, '*488nm*.tif*'));
list = {list.name};

%% iterate through the files
% ignore warnings for unknown tags for current session
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

nFile = numel(list);
tOuter = tic;
for iFile = 1:nFile
    tInner = tic;
    
    %% load the file
    fileName = list{iFile};
    fprintf('[%s]\n', fileName);
    filePath = fullfile(srcDir, fileName);
    
    I = tiff.TIFFStack(filePath);
    [nx, ny, nz] = size(I);
    I = I(:);
    I = reshape(I, nx, ny, nz);
    
    %% pre-process
    % re-order the stack to phase-wise
    
    
    if iFile == 1
        %% find out kp value using the first stack
    end
    
    %% reconstruction
    
    %% cleanup and statistics
    fclose('all');
    
    tElapse = toc(tInner);
    fprintf('%f seconds elapsed\n', tElapse);
end
tElapse = toc(tOuter);
fprintf('%f seconds to process the incoming data\n', tElapse);
