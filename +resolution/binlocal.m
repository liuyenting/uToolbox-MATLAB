function [I, pxsz] = binlocal(coords, sz)
%BINLOCAL Bin the localization data
%
%   TBA

% estimate proper pixel dimensions that can contain all the data
pxsz = estpxsize(coords, sz);

% generate the empty stack
I = zeros(sz);

% convert coordinates to pixel scale
pxsz = repmat(pxsz, [size(coords, 1), 1]);
coords = coords ./ pxsz;
coords = round(coords) + 1;

% fill them to the stack
ind = ndsub2ind(sz, coords);
for i = ind
    I(i) = I(i)+1;
end

end

function pxsz = estpxsize(coords, sz)
%ESTPXSIZE Estimate the pixel size of the binned image.
%
%   PXSZ = ESTPXSIZE(COORDS, SZ) estimates the generated pixel size 
%   according to the maxmium size of the coordinate range.

maxcoord = max(coords);
pxsz = maxcoord ./ (sz-1);
pxsz = repmat(max(pxsz(:)), size(sz));

end
