clear all; close all; %#ok<CLALL>

sz = [1024, 1024];
r = 10;

I = zeros(sz);

% find the midpoint
midpt = sz / 2;

pres = 1;
% generate the sampling location
ares = calcares(pres, r);
angles = 0:ares:2*pi;

% generate the Cartesian coordinates
[x, y] = pol2cart(angles, r);

% generate comparison grid
[xi, yi] = meshgrid(1:sz(1), 1:sz(2));
xi = xi - midpt(1);
yi = yi - midpt(2);

ri = rdist(xi, yi);
I((ri >= r) & (ri < r+1)) = 1;

figure;
imagesc(I);
    axis image;

function r = rdist(x, y)

r = sqrt(x.^2 + y.^2);

end

function ares = calcares(pres, r)
%CALCARES Calculate angular resolution from desired pixel resolution.

ares = pres / r;

end