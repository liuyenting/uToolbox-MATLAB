function [npx, newpxsize] = estsize(coords, pxsize, mag)
%ESTSIZE Estimate super-resolved image size.

% minimal original pixel to describe the coordinates
maxcoord = max(coords);
npx = ceil(maxcoord ./ pxsize);

% magnifiy
npx = npx .* mag;

% new pixel size
newpxsize = pxsize ./ mag;

end
