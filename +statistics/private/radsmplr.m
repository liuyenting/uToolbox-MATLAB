function A = radsmplr(ip, midpt, r, pres)
%RADSMPLR Sampling data in radial direction.
%
%   A = RADSMPLR(IP, MIDPT, R, PRES) samples at radius R with pixel
%   resolution PRES using IP and MIDPT as the center. Output array A is the
%   sampled data.

% angular resolution from desired pixel resolution
ares = pres/r;
% sampling locations
angles = 0:ares:2*pi;

% generate the Cartesian coordinates
[x, y] = pol2cart(angles, r);
% sample on the matrix
A = ip(x+midpt(1), y+midpt(2));

end
