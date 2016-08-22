clearvars;
%% Parameters
NA = 1.44;
lambda = 488;   % nm
ni = 1.3;

axial = (lambda/8) * 1e-3;  % um
radial = 0.01;   % um

nx = 96;        % pixels
ny = 96;        % pixels
nz = 87;       % pixels

%% Create the PSF.
psf = BornWolfPSF();

psf.setopticsparam(NA, lambda, ni);
psf.setresparam(axial, radial);
psf.setoutputparam([nx, ny, nz]);

psf.generate();

%% Plot the result.
%close all;

psf.plotortho('xz');
