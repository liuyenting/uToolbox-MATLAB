function TF = psf2tf(imSz, PSF, M, kp, parms)
%PSF2TF Convert PSF to transfer functions.
%   
%   TBA

%% parameters
nOri = parms.Orientations;
nPhase = parms.Phases;

% identify the type of SIM
%   * 2-D SIM requires 3 phases
%   * 3-D SIM requires 5 phases
if nPhase == 3
    nd = 2;
elseif nPhase == 5
    nd = 3;
else
    error(generatemsgid('InvalidInput'), ...
          'Unsupported type of SIM data source.');
end

psfSz = size(PSF);

pxSz = parms.PixelSize;

n = parms.RefractiveIndex;
lambda = parms.Wavelength;

%% pre-allocate
% domains
D = zeros([imSz, nPhase], 'single');
% the transfer functions
TF = zeros([imSz, nPhase], 'single');

%% generate illumination patterns
% use the last two components for frequency estimation (kx, kz)
%   * 2-D SIM m1
%   * 3-D SIM m2
kp = kp(:, end-1:end, :);

% convert from pixel to spatial frequency
kp = kp ./ (imSz.*pxSz).';
% combine frequency components
f = hypot(kp(1, :, :), kp(2, :, :));
% average the frequency values from each orientation and -/+ terms
f = reshape(f, 2*nOri, 1);
f = mean(f);

% turns into the actual wave number
%   k = 2*pi / p
%     = 2*pi * f
% Note: While p is period and f is frequency, we are working backward from 
% the reciprocal space, hence f is used instead of p.

% beam angle
%   p = lambda / (2*n*sin(theta)), p is period
%   sin(theta) = lambda / (2*n*p)
%              = lambda / (2*n) * f

% kx*sin(theta) = (2*pi * f) * lambda / (2*n) * f
%               = 2*pi * f^2 * lambda / (2*n)
%               = A
A = 2*pi * f^2 * lambda / (2*n);

% phase array
phi = linspace(0, 2*pi, nPhase+1);
phi = phi(1:end-1).';
% insert singleton dimension to enable implicit padding
phi = reshape(phi, 1, 1, []);

% the grid for computation
[vx, ~] = meshgrid(1:psfSz(1), 1:psfSz(2));
% convert to real word scale
vx = vx*pxSz(1);

if nd == 2
    % I = 2*cos(2*kx*sin(theta) + phi)
    Im = 2*cos(2*A*vx + phi);
elseif nd == 3
    % kz component
    %   fz = (1-cos(theta)) * n/lambda
    %      = (1-sqrt(1-sin(theta)^2)) * n/lambda
    %      = n/lambda - sqrt((n/lambda)^2 - f^2)
    %   kz = 2*pi * fz
    kz = 2*pi * (n/lambda - sqrt((n/lambda)^2 - (f/2)^2));

    % I = 4*cos(kx*sin(theta) + phi)*cos(kz*(cos(theta)-1)) + 
    %     2*cos(2*kx*sin(theta) + 2*phi) 
    %
    % kz*(cos(theta)-1) = kz*(-kz * lambda/n)
    %                   = - kz^2 * lambda/n
    %                   = B
    B = - kz^2 * lambda/n;
    
    % Since we are working with planes instead of a volume, z is always 1
    % for the kz term (B).
    Im = 4*cos(A*vx + phi)*cos(B) + 2*cos(2*A*vx + 2*phi);
end
% shift lowest value from negative to 0
% Note: 2-D SIM has 2 beam sources, while 3-D SIM has 3.
Im = Im + nd;

if parms.Debug
    figure('Name', 'Illumination Patterns', 'NumberTitle', 'off');
    for iPhase = 1:nPhase
        subplot(1, nPhase, iPhase);
        imagesc(Im(:, :, iPhase));
            axis image;
            title(sprintf('Phase %d', iPhase));
    end
end

%% retrieve domains (bands separation)
% apply the modulations
PSF = PSF .* Im;

% retrieve the reciprocal space images
for iPhase = 1:nPhase
    D(:, :, iPhase) = psf2otf(PSF(:, :, iPhase), imSz);
end

% flatten the array
D = reshape(D, [prod(imSz), nPhase]);
% solve the matrix
D = (M \ D')';
% reshape back to original image size
D = reshape(D, [imSz, nPhase]);

if parms.Debug
    figure('Name', 'Transfer Functions', 'NumberTitle', 'off');
    for iPhase = 1:nPhase
        subplot(1, nPhase, iPhase);
        imagesc(abs(ifftshift(D(:, :, iPhase))));
            axis image;
    end
end

%% remove initial phase k0
for iPhase = 2:2:nPhase
    Dm = D(:, :, iPhase);
    Dp = D(:, :, iPhase+1);
    
    options = optimoptions(@lsqnonlin, ...
                           'Display', 'iter-detailed', ...
                           'FiniteDifferenceType', 'central', ...
                           'FiniteDifferenceStepSize', 2*pi/imSz(1));
    s = lsqnonlin(@(s)(errfunc(s, Dm, Dp)), pi, [], [], options);
    
    D(:, :, iPhase) = exp(-1i*s)*Dm;
    D(:, :, iPhase+1) = exp(1i*s)*Dp;
end

    % cost function for the non-linear fitting
    function err = errfunc(s, m, p)
        err = exp(-1i*s)*m - exp(1i*s)*p;
        % return type is required to be double
        err = double(err);
    end

%% bead mask
%TODO compensate for non-deal point source

%% average and LPF
% cutoff frequency
f = 2*parms.NA / lambda;
% cut-off radius in pixel
r = f * (imSz.*pxSz);

midpt = floor(imSz/2)+1;
for iPhase = 1:nPhase
    T = D(:, :, iPhase);
    % shift to natural order for processing
    T = ifftshift(T);
    
    % create radial profile
    % Note: Due to the limitation of radial sampler, only square image is
    % plausible for now.
    TF(:, :, iPhase) = statistics.radialmean(T, midpt, floor(r));
end

end
