function sz = estimsize(coords, res, varargin)
%ESTIMSIZE Estimate the size of the binned image.
%
%   SZ = ESTIMSIZE(COORDS, RES) estimates the generated pixel size 
%   according to the maxmium size of the coordinate range.

maxcoord = max(coords);
sz = maxcoord / res;
sz = ceil(sz);

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
