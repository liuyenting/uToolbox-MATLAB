function s = radialsum(I, r, res)
%RADIALSUM Calculate the radial sum of an image with interpolation.
%
%   S = RADIALSUM(I, R)
%   S = RADIALSUM(I, R, RES) calculate radial sum at R spatial frequency
%   using resolution of RES. RES is default to 1 degree if not assigned.

% ensure the input is correct
if ~ismatrix(I)
    error('resolution:radialsum', 'Not a 2-D data.');
end

if nargin == 2
    res = 1;
end

% find the midpoint
midpoint = size(I) / 2;

% generate the sampling location
angles = 0:res*pi/180:pi;
angles = angles';
% convert to Cartesian coordinates with specified frequency scale
pq = r * [cos(angles), sin(angles)];
% offset
pq = pq + repmat(midpoint, [size(pq, 1), 1]);

s = interp2(I, pq(1), pq(2));
s = sum(s);

end
