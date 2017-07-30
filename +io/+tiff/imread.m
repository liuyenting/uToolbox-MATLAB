function varargout = imread(path, nowarn)

%% Suppress warnings from libtiff.
if nargin == 1
    nowarn = false;
end
if nowarn
    % Ignore warnings for unknown tags.
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
end

%% Load the stack.
imgStack = tiff.TIFFStack(path);
[nx, ny, nz] = size(imgStack);
if nz == 1
    % use standard imread to bypass bug in TIFFStack
    clear imgStack;
    img = imread(path);
else
    img = imgStack(:);
    img = reshape(img, nx, ny, nz);
end
varargout{1} = img;

%% Save the stack size if needed.
if nargout == 2
    varargout{2} = [nx, ny, nz];
elseif nargout == 4
    varargout{2} = nx;
    varargout{3} = ny;
    varargout{4} = nz;
elseif nargout ~= 1
    warning(generatemsgid('TooManyOutput'), ...
            'Invalid output arguments schema, only data is returned.');
end

%% Re-enable the libtiff warnings.
if nowarn
    warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
end

end