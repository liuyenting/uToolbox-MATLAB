function I = binlocal(coords, npx, pxsz)
%BINLOCAL Bin the localization data
%
%   TBA

% generate the empty stack
I = zeros(npx);

% convert coordinates to pixel scale
pxsz = repmat(pxsz, [size(coords, 1), 1]);
coords = coords ./ pxsz;
coords = round(coords);

% fill them to the stack
ind = ndsub2ind(npx, coords);
I(ind) = 1;

end

function ind = ndsub2ind(sz, sub)
%NDSUB2IND Convert N dimensional subscripts to linear indices.

% cumulative products for the dimension
k = cumprod(sz);
k = [1, k(1:end-1)];

ind = sub * k';

end
