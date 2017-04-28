function [npx, newpxsize] = estsize(coords, pxsize, mag)
%ESTSIZE Estimate super-resolved image size.
%
%   [NPX, NEWPXSIZE] = ESTSIZE(COORDS, PXSIZE, MAG) TBA
%
%   Note: This function is based on the original pixel size, so the
%   estimated result might be greater than the smallest size it can be.

% minimal pixel to describe the coordinates
maxcoord = max(coords);
npx = ceil(maxcoord ./ pxsize);

% magnifiy
npx = npx .* mag;
% additional spacing
npx = npx + 1;

% new pixel size
newpxsize = pxsize ./ mag;

end
