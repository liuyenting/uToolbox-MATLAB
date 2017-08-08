function sheardemo(file, acqParam)
%SHEARDEMO Demonstrate the basic usage of class Shear.
%   SHEARDEMO
%   SHEARDEMO(FILE) will shear the input stack according to presets. If
%   FILE is not assigned, the script will ask for one.
%
%   SHEARDEMO(FILE, ACQPARAM) will shear the input stack according to
%   ACQPARAM instead of presets.

%% Parameters.
if nargin == 0
    [oldImgName, imgDir] = uigetfile('*.tif*', ...
                                     'Pick a file', ...
                                     'MultiSelect', 'off');
    if imgDir == 0
        % Cancelled.
        return;
    else
        oldImgPath = fullfile(imgDir, oldImgName);
    end
elseif nargin > 0
    oldImgPath = file;
end
    
if nargin < 2
    % Acquisition paramters.
    acqParam.ObjectiveAngle = 32.8; % [deg]
    acqParam.ZStepWidth = 0.5;       % [um]
    acqParam.PixelWidth = 0.102;     % [um]
end

%% Read TIFF stack into host memory.
tic;

% Ignore warnings for unknown tags.
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

oldImgStack = tiff.TIFFStack(oldImgPath);
[nx, ny, nz] = size(oldImgStack);
% Preallocate host memory.
oldImg = oldImgStack(:);
oldImg = reshape(oldImg, nx, ny, nz);

warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

tElapse = toc;
fprintf(' [[[ %f seconds to load the raw image ]]]\n', tElapse);

%% Begin processing.
tic;

shObj = postproc.Shear();
shObj.setacqparam(acqParam);

shObj.loadstack(oldImg);
shObj.execute();
newImg = shObj.retrieveresult();

tElapse = toc;
fprintf(' [[[ %f seconds to complete the shearing ]]]\n', tElapse);

%% Save the result.
tic;

newImg = uint16(newImg);
[~, fileName, fileExt] = fileparts(oldImgName);
newImgPath = [imgDir, fileName, '_sheared', fileExt]; 
io.fwrite(newImg, newImgPath, true);

tElapse = toc;
fprintf(' [[[ %f seconds to save the result ]]]\n', tElapse);

end