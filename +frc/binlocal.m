function I = binlocal(coord, pxsz, mag)
%BINLOCAL Bin the localization data
%
%   BUILDSUPER(COORD, PXSZ, MAG) constructs a binned stack image of
%   localization coordinates COORD, scaled by the original camera pixelsize
%   PXSZ in nm, with an magnification ratio MAG. Typically, MAG is
%   greater than 4.

coord = offsetorigin(coord);

% calculate the super-resolved image size in pixel
npx = pxrange(coord, pxsz);
npx = npx .* mag;

% generate the empty stack
I = zeros(npx);

% convert pixel size to super-resolved scale
pxsz = pxsz / mag;
% convert coordinates to pixel scale
coord = coord ./ pxsz;
coord = round(coord);

min(coord(:, :))
max(coord(:, :))

% fill them to the stack
ind = ndsub2ind(npx, coord);
I(ind) = 1;

end

function newcoord = offsetorigin(coord)
%OFFSETORIGIN Move the coordinates to starts from 0-based origin.

mincoord = min(coord);
n = size(coord, 1);
mincoord = repmat(mincoord, [n, 1]);
newcoord = coord - mincoord;

end

function npx = pxrange(coord, pxsz)
%PXRANGE Calculate the minimal pixels to describe the coordinate range.

maxcoord = max(coord);
npx = ceil(maxcoord ./ pxsz);

end

function ind = ndsub2ind(sz, sub)
%NDSUB2IND Convert N dimensional subscripts to linear indices.

% cumulative products for the dimension
k = cumprod(sz);
k = [1, k(1:end-1)];

ind = sub * k';

end
