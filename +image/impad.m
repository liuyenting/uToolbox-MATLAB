function [J, varargout] = impad(I, sz, varargin)
%IMPAD Pad the image.
%
%   J = IMPAD(I, SZ) pad the image to specified size, 0s are padded around
%   the image.
%   J = IMPAD(I, SZ, METHOD) controls the location of the padded image,
%   method can be
%       'center'    Image locates in the center.
%       'corner'    Image lcoates at the lower-left corner in matrix-wise.
%   Default method is 'corner'.
%   J = IMPAD(..., 'PadValue', V) pads the surrounding region use value V,
%   V is default to 0.
%   [J, C] = IMPAD(...) returns the anchor of the original image in C.

imsz = size(I);
% ignore the request if sizes are matched
if imsz == sz
    J = I;
    if nargout == 2
        varargout{1} = [1, 1];
    end
    return;
end

% size differences.
dsz = sz - imsz;
if any(dsz < 0) 
    error('image:impad', ...
          'Output size has to be greater or equal to the input size.');
end

% parse the additional parameters
p = inputParser;
addOptional(p, 'Method', 'corner', @(x)(ischar(x)));
addParameter(p, 'PadValue', 0);
parse(p, varargin{:});

method = p.Results.Method;
padval = p.Results.PadValue;

% create the output image
if padval ~= 0
    J = padval*ones(sz);
else
    J = zeros(sz);
end

% calculate the corner
if strcmp(method, 'corner')
    c = [0, 0];
elseif strcmp(method, 'center')
    c = floor(dsz/2);
else 
    error('image:impad', 'Unknown padding method.');
end
c = c + [1, 1];

% paste the array
J(c(1):c(1)+imsz(1)-1, c(2):c(2)+imsz(2)-1) = I;

if nargout == 2
    varargout{1} = c;
end

end