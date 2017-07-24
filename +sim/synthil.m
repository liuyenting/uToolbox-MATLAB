function IL = synthil(imSz, kp, parms)
%SYNTHIL Synthesize the illumination pattern.
%   Detailed explanation goes here

%% parameters
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
f = mean(f(:));
fprintf('\tSynthesized Pattern Period = %.4fnm\n', 1/f);

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
    IL = 2*cos(2*A*vx + phi);
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
    IL = 4*cos(A*vx + phi)*cos(B) + 2*cos(2*A*vx + 2*phi);
end
% shift lowest value from negative to 0
% Note: 2-D SIM has 2 beam sources, while 3-D SIM has 3.
IL = IL + nd;

if parms.Debug
    figure('Name', 'Illumination Patterns', 'NumberTitle', 'off');
    for iPhase = 1:nPhase
        subplot(1, nPhase, iPhase);
        imagesc(IL(:, :, iPhase));
            axis image;
            title(sprintf('Phase %d', iPhase));
    end
end

end

