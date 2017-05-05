function [I, pxsz] = binlocal(coords, sz, pxsz)
%BINLOCAL Bin the localization data
%
%   TBA

if size(coords, 2) ~= length(sz)
    error('resolution:binlocal', ...
          'Dimension of the image size mismatch the coordinate list.');
end

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
