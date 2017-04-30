function [I, pxsz] = binlocal(coords, sz, pxsz)
%BINLOCAL Bin the localization data
%
%   TBA

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
