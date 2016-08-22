function localizedemo
%LOCALDEMO Demonstrate the basic usage of class Localize.
%   LOCALDEMO will try to perform localization on presented image stack.

%% Presets.
imgFilename = 'local_spl1.tif';
imgDir = 'data';

%% Retrieve the image.
imgPath = fullfile(imgDir, imgFilename);

tic;

% Ignore warnings for unknown tags.
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

imgStack = tiff.TIFFStack(imgPath);
[nx, ny, nz] = size(imgStack);
% Preallocate host memory.
imgData = imgStack(:);
imgData = reshape(imgData, nx, ny, nz);

warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

tElapse = toc;
fprintf(' [[[ %f seconds to load the raw image ]]]\n', tElapse);

end