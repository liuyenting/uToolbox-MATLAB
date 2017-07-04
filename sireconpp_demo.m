%SIRECONPP_DEMO Demonstrate the plan-by-plan reconstruction feature.
%
%   TBA

close all;
clearvars;

%TODO remove OS dependent test code
if ispc
    srcDir = 'F:\\Andy\\06302017_SIM\\Tan1_3_3DSIM';
    psfFileName = 'F:\\Andy\\06302017_SIM\\PSF\\PSF_NA1p1_520nm_x100nm.tif';
else
    srcDir = 'data/sim/06302017_SIM/Tan1_3_3DSIM';
    psfFileName = 'data/sim/06302017_SIM/PSF/PSF_NA1p1_520nm_x100nm.tif';

    % prepend the root
    srcDir = fullfile(userpath, srcDir);
    psfFileName = fullfile(userpath, psfFileName);
end

dstSuffix = '_sirecon';
dstDir = [srcDir, dstSuffix];

% orientations
nOri = 1;
% phases
nPhase = 5;

%% verify the input
% check whether the input directory exsists
if exist(srcDir, 'dir') ~= 7
    error('sim:sireconpp_demo', 'Input directory does not exist.');
end

% check whether the output directory is empty
if exist(dstDir, 'dir') == 7
    content = dir(dstDir);
    if numel(content) > 2
        error('sim:sireconpp_demo', 'Output directory is not empty.');
    end
else
    % create the directory
    status = mkdir(dstDir);
    if ~status
        error('sim:sireconpp_demo', ...
              'Unable to create the output directory.');
    end
end

%% load the file list
%TODO only load 488nm for simplicity
list = dir(fullfile(srcDir, '*488nm*.tif*'));
list = {list.name};

%% preload the PSF
Ipsf = tiff.imread(psfFileName);
Ipsf = single(Ipsf);

% center the PSF
Ipsf = image.centerpsf(Ipsf);

% normalize PSF, ensure kernel sum is 1
Ipsf = Ipsf / sum(Ipsf(:));

%% iterate through the files
% ignore warnings for unknown tags for current session
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

nFile = numel(list);
tOuter = tic;
%OVERRIDE
nFile = 1;
for iFile = 1:nFile
    tInner = tic;
    
    %% load the file
    fileName = list{iFile};
    fprintf('[%s]\n', fileName);
    filePath = fullfile(srcDir, fileName);
    
    I = tiff.TIFFStack(filePath);
    sz = size(I);
    
    %% pre-process
    % re-order the stack to phase-wise
    I = sim.phmajor(I, sz, nOri, nPhase);
    
    if iFile == 1
        %% create calibration preview
        calDir = fullfile(dstDir, 'calibrate');
        
        % create the directory
        status = mkdir(calDir);
        if ~status
            error('sim:sireconpp_demo', ...
                  'Unable to create the calibration output directory.');
        end

        % create MIP for debug
        Ip = zeros([nOri, nPhase, sz(1:2)]);
        for iOri = 1:nOri
            for iPhase = 1:nPhase
                % extract the volume
                P = I(iOri, iPhase, :, :, :);
                P = squeeze(P);
                
                % MIP along Z axis
                P = max(P, [], 1);
                P = squeeze(P);
                
                % save the output
                fname = sprintf('o%d_p%d.tif', iOri, iPhase);
                fpath = fullfile(calDir, fname);
                tiff.imsave(P, fpath);
                
                % save the result for further processing
                Ip(iOri, iPhase, :, :) = P;
            end
        end
        
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
