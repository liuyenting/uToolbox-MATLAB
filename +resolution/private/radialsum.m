function s = radialsum(I, smplratio)
%RADIALSUM Calculate the radial sum of an image with interpolation.
%
%   S = RADIALSUM(I, ARES) calculate radial sum at different spatial 
%   frequency using angular resolution of ARES. ARES is default to 1 degree
%   if not assigned.
%   SMPLRATIO implies the sampling ratio, (0, 1], default at 0.5.

% find the minimal dimension
sz = size(I);
if sz(1) ~= sz(2)
    warning('resolution:radialsum', ...
            'Not a square matrix, cropped to the minimal square.');
    [I, sz] = sqcrop(I);
end

if nargin == 1
    % default cutoff at 50% (Nyquist criterion)
    smplratio = 0.5;
else
    if (smplratio <= 0) || (smplratio > 1)
        error('resolution:radialsum', ...
              'Sampling ratio can only range from (0, 1].');
    end
end

% radius sample location
r = 1:floor((sz(1)/2) * smplratio);
% find the center
[~, ind] = max(I(:));
[xi, yi] = ind2sub(sz, ind);
% generate the distance map
rm = distmap(sz, [xi, yi]);

% start sampling
s = zeros(size(r));
for i = 1:length(s)
    ind = (rm >= r(i)) & (rm < r(i)+1);
    s(i) = sum(I(ind));
end

end

function r = distmap(sz, midpt)
%DISTMAP Generate distance map from the center of the image.

if nargin == 1
    % find the midpoint
    midpt = sz / 2;
end

% generate comparison grid
[xq, yq] = meshgrid(1:sz(1), 1:sz(2));
xq = xq - midpt(1);
yq = yq - midpt(2);

r = sqrt(xq.^2 + yq.^2);

end

function [Iout, sz] = sqcrop(Iin)
%SQCROP Crop input array to square array.

[nrow, ncol] = size(Iin);
if nrow == ncol
    Iout = Iin;
    sz = [nrow, ncol];
else
    sz = min(nrow, ncol);
    redsz = max(nrow, ncol) - sz;
    
    % padded size is half of the redundant size
    padsz = floor(redsz/2);
    
    if sz == nrow
        Iout = Iin(:, padsz:padsz+sz-1);
    else
        Iout = Iin(padsz:padsz+sz-1, :);
    end
end

end