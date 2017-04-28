function [A, lbl] = radialsmplr(I, rstep, ares)
%RADIALSMPLR Sampling data in radial direction.
%
%   TBA

% ensure the input is correct
if ~ismatrix(I)
    error('resolution:radialsum', 'Not a 2-D data.');
end

L = size(I, 1);
if size(I, 2) ~= L
    warning('resolution:radialsum', 'Not a square matrix.');
end

% find the midpoint
midpt = size(I) / 2;

% generate the sampling location
angles = 0:ares*pi/180:pi;
angles = angles';

% generate the Cartesian coordinates
pq = rstep/L * [cos(angles), sin(angles)];
% shift the center
pq = pq + repmat(midpt, [size(pq, 1), 1]);

% sample on the matrix
A = interp2(I, pq(1), pq(2));

if nargout == 2
    lbl = genlbl(pq, size(I));
end

end

function lbl = genlbl(pts, npx)
%GENLBL Generate the label file.

lbl = zeros(npx);

% fill the points on the label matrix
pts = round(pts);
ind = ndsub2ind(npx, pts);
lbl(ind) = 1;

end

