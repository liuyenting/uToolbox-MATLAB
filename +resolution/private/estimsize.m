function sz = estimsize(coords, res, varargin)
%ESTIMSIZE Estimate the size of the binned image.
%
%   SZ = ESTIMSIZE(COORDS, RES) estimates the generated binned image size
%   according to the specified resolution.
%   SZ = ESTIMSIZE(..., PARAM) allows one to specify croping or extending
%   method.
%
%   Parameters
%   ----------
%   'MaxSquare'     Extend the region to a square that can contain the 
%                   entire image. Rest of the region is filled with zero.
%   'MinSquare'     Crop the region to a square with its length the minimum
%                   of the rectangle. 

maxcoord = max(coords);
sz = maxcoord / res;
% unconditional carry
sz = ceil(sz + 1);

if nargin == 3
    % additional resolution setting
    if nargin > 3
        warning('resolution:estimsize', ...
                'Excessive arguments are ignored.');
    end
    type = varargin{1};

    if strcmp(type, 'MaxSquare')
        sz = repmat(max(sz(:)), size(sz));
    elseif strcmp(type, 'MinSquare')
        sz = repmat(min(sz(:)), size(sz));
    else
        error('resolution:estimsize', 'Unknown resolution arguments.');
    end
end

end
