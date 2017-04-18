function np = npcgen(f, ori, r, unc, ccd, cutoff, psf, snr)
% NPCGEN Generates simulated NPC dSTORM camera view.
%
%   NP = NPCGEN(F, ORI, R, UNC, CCD, CUTOFF, PSF, SNR) 
%       F       Figure handles, contains 'Preview' and 'Camera'
%       ORI     [azimuth, elevation] 1x2 array, [deg]
%       R       Radius of the NPC octagon, [nm]
%       UNC     [radial uncertainty, axial uncertainty], [nm]
%       CCD     CCD configuration, 'XDim' and 'YDim' are in pixels,
%               'PixelSize' in [nm].
%       CUTOFF  Control the blinking probability, (0, 1).
%       PSF     Describe an ellipsoid, 'RadialSigma' and 'AxialSigma' in
%               [nm], 'MaxInt' is the maximum intensity in [au]. 'SmplNSig'
%               control the sampling size. 
%       SNR     Background noise SNR, in linear scale.

%% extract parameters
if isempty(f)
    nofig = 1;
else
    nofig = 0;
    h_pre = f.Preview;
    h_cam = f.Camera;
end

% azimuth [deg]
az = ori(1);
% elevation [deg]
el = ori(2);

% ccd dimension, (x, y) [px]
ccdDim = [ccd.XDim, ccd.YDim];
% pixel size [nm]
px = ccd.PixelSize;

% fluorophore intensity [au]
maxInt = psf.MaxIntensity;
% radial sigma [nm]
sig_xy = psf.RadialSigma;
% axial sigma [nm]
sig_z = psf.AxialSigma;
% psf size, n sigma [#]
nsig = psf.SmplNSig;

%% flags
% convolve psf?

% background noise?

%% generate the octagon
% default position
p = genoct(r);

% rotate elevation
p = rotx(el) * p;
% rotae azimuth
p = rotz(az) * p;

% blink the fluorophores
bp = blink(p, cutoff);
% apply uncertainty
if ~isempty(unc)
    unc_xy = unc(1);
    unc_z = unc(2);
    up = uncert(bp, unc_xy, unc_z);
else
    up = bp;
end

if ~nofig
    figure(h_pre);

    % plot the dots
    subplot(1, 2, 1);
    % plot absolute positions
    scatter3(p(1, :), p(2, :), p(3, :), 'filled', 'MarkerFaceColor', 'b');
    axis([-r, r, -r, r, -r, r]);
    axis square;
    xlabel('x', 'FontWeight', 'bold'); 
    ylabel('y', 'FontWeight', 'bold'); 
    zlabel('z', 'FontWeight', 'bold');
    view(60, 30);
    % plot pertubed positions
    hold on;
    scatter3(up(1, :), up(2, :), up(3, :), 'filled', 'MarkerFaceColor', 'r');
    % annotate the viewing angle
    str = ['Az = ' num2str(az) '^{\circ}, El = ' num2str(el) '^{\circ}'];
    title(str);
    % plot projection direction
    quiver3(0, 0, 0, 0, 0, -r, 'MaxHeadSize', 1, 'Color', 'k');
    hold off;

    % plot camera projection
    subplot(1, 2, 2);
    % plot absolute positions
    scatter(p(1, :), p(2, :), 'filled', 'MarkerFaceColor', 'b');
    axis([-r, r, -r, r]);
    axis square;
    xlabel('x', 'FontWeight', 'bold'); 
    ylabel('y', 'FontWeight', 'bold'); 
    title('Projection (0, 0, -1)');
    % plot pertubed positions
    hold on;
    scatter(up(1, :), up(2, :), 'filled', 'MarkerFaceColor', 'r');
    hold off;
end

%% convert to camera grid
% interpolate to the camera grid
vq = intpcam(ccdDim, px, maxInt, sig_z, up);

% convolve psf
psf = psfgen(px, sig_xy, nsig);
vq = conv2(vq, psf, 'same');

% add thermal noise
np = noisy(vq, snr);

if ~nofig
    figure(h_cam);
    imagesc(np);
    set(gca, 'YDir', 'normal');
    axis image;
    colormap(gray);

    drawnow;
end

end

function p = genoct(r)
% GENOCT Generate an octagon of specific radius.
%
%   P = GENOCT(R) returns a point array P of an octagon with their radius
%   set to R.

a = linspace(0, 2*pi, 9);
a = a(1:8);

[x, y] = pol2cart(a, r);
z = zeros(size(x));

p = [x; y; z];

end

function m = rotx(d)
% ROTX Rotation matrix for rotations around x-axis.

r = deg2rad(d);

m = [
    1, 0, 0;
    0, cos(r), -sin(r);
    0, sin(r), cos(r)
];

end

function m = roty(d)
% ROTY Rotation matrix for rotations around y-axis.

r = deg2rad(d);

m = [
    cos(r), 0, sin(r);
    0, 1, 0;
    -sin(r), 0, cos(r)
];

end

function m = rotz(d)
% ROTZ Rotation matrix for rotations around z-axis.

r = deg2rad(d);

m = [
    cos(r), -sin(r), 0;
    sin(r), cos(r), 0;
    0, 0, 1
];

end

function p = blink(p, cutoff)
% BLINK

np = size(p, 2);

% generate random level
rng('shuffle');
lv = rand(1, np);

% extract the index
ind = lv < cutoff;
% remove columns by index
p(:, ind) = [];

end

function po = uncert(pi, uxy, uz)
% NOISY Apply uncertainties to input coordiantes.
%
%   PO = UNCERT(PI, UXY, UZ) applies uncertainty UXY on X and Y axis, while
%   uncertinaty UZ on Z axis to data points PI. The fluctuation range is
%   (-UXY, UXY) and (-UZ, UZ) respectively.

% pre-allocation
po = zeros(size(pi));
np = size(pi, 2);

% perturbe x and y
rng('shuffle');
po(1:2, :) = pi(1:2, :) + (2*uxy * rand(2, np) - uxy);
% perturbe z
rng('shuffle');
po(3, :) = pi(3, :) + (2*uz * rand(1, np) - uz);

end

function [vx, vy] = gencam(sz, px)
% GENCAM Generate camera sensor area.

rx = genrwvec(sz(1), px);
ry = genrwvec(sz(2), px);
[vx, vy] = meshgrid(rx, ry);

end

function v = genrwvec(sz, px)
% GENRWVEC Generate real world scale vector set.

% lower and upper bound in pixel
ub = (sz-1)/2;
lb = -ub;

v = linspace(lb, ub, sz) * px;

end

function vq = intpcam(dim, px, mi, sig, p)
% INTPCAM Interpolates scattered data to camera grid evenly.
%
%   VQ = INTPCAM(DIM, PX, MI, SIG, P) spread data points P onto a CCD of 
%   dimension DIM with pixel size PX as VQ. MI is the maximum instensity of
%   the fluorophore, while SIG is the sigma of the PSF on axial direction.
%
%   Note: Using linear interpolation.

offset = (dim-1)/2 * px;
% result image
vq = zeros(dim);

% shift z, zero at the shallow one
z = p(3, :);
z = z - max(z);
int = zint(mi, sig, z);

np = size(p, 2);
for i = 1:np
    % real world coordinate
    rx = p(1, i) + offset(1);
    ry = p(2, i) + offset(2);
    fx = rx/px;
    fy = ry/px;
    ix = floor(fx);
    iy = floor(fy);
    
    % disperse them to neighbors
    vq(iy:iy+1, ix:ix+1) = disperse(int(i), fx, fy);
end

end

function i = zint(mi, sig, z)
% ZINT Attenuate max intensity by depth.
%
%   I = ZINT(MI, SIG, Z) TBA

i = mi * exp(-(z/sig).^2 / 2);

end

function dm = disperse(v, x, y)
% DISPERSE Disperse non integer coordinate to its neighbor pixels.
%
%   DM = DISPERSE(V, X, Y) disperses value V to non-integer coordinate 
%   (X, Y)'s neighbors. DM is a 2x2 matrix that describes the neighbor
%   pixels.

% pre-allocation
dm = zeros(2, 2);

x0 = floor(x);
x1 = ceil(x);

y0 = floor(y);
y1 = ceil(y);

% NOTE: fraction values are crossed over, since the close the position is,
% the greater the weighting should be.
if x0 == x1
    % x singleton
    if y0 == y1
        % direct assignment
        dm(1, 1) = v;
    else
        % disperse y
        vu = v * (y-y0);
        vd = v * (y1-y);
    
        dm(1:2, 1) = [vd; vu];
    end
else
    % left-right (x)
    vl = v * (x1-x);
    vr = v * (x-x0);

    if y0 == y1
        % y singleton
        dm(1, 1:2) = [vl, vr];
    else
        % full disperse
        vlu = vl * (y-y0);
        vld = vl * (y1-y);
        vru = vr * (y-y0);
        vrd = vr * (y1-y);

        dm(1:2, 1:2) = [
            vld, vrd;
            vlu, vru
        ];
    end
end

end

function psf = psfgen(px, sig, nsig)
% PSFGEN Generate a simple Gaussian point spread function.

% determine size of the grid
npx = 2 * sig/px;
vt = linspace(-nsig*sig, nsig*sig, npx);
% generate the grid itself
[vx, vy] = meshgrid(vt);

% generate the psf function
psf = exp(-(vx.^2 + vy.^2)/(2 * sig^2));

end

function po = noisy(pi, snr)
% NOISY Apply white noise to the background.
%
%   PO = NOISY(PI, SNR) apply a background noise of signal-to-noise ratio
%   SNR to the input image PI.

% probe for size and max of the input
sz = size(pi);
mi = max(pi(:));
% calculate the noise level
nl = mi / snr;

% generate the noise profile
no = nl .* rand(sz);
% apply to the data
po = pi + no;

end
