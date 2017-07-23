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
hIL = figure('Name', 'Illumination Pattern', 'NumberTitle', 'off');
hPost = figure('Name', 'Transfer Function', 'NumberTitle', 'off');

% domains
D = zeros([imSz, nPhase], 'single');
% the transfer functions
TF = zeros([imSz, nPhase], 'single');

% the grid for computation
[vx, ~] = meshgrid(1:psfSz(1), 1:psfSz(2));
% convert to real word scale
vx = vx * pxSz(1);

%% generate illumination patterns
% use the last two components for frequency estimation (kx, kz)
%   * 2-D SIM m1
%   * 3-D SIM m2
kp = kp(:, end-1:end, :);

% convert units of kp values to spatial frequency
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

if nd == 2
    % I = 2*cos(2*kx*sin(theta) + phi)
    %TODO
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
    figure(hIL);
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

%DEBUG
figure;
for iPhase = 1:nPhase
    subplot(1, nPhase, iPhase);
    imagesc(abs(ifftshift(D(:, :, iPhase))));
        axis image;
end

%% process
for iOri = 1:nOri
    for iPhase = 2:nPhase
        %% calculate the modulation frequency by kp
        k = kp(:, iPhase-1, iOri);
        % convert to spatial frequency
        k = k ./ (imSz.*pxSz).';

        % frequency of the diffraction pattern
        f = hypot(k(1), k(2));
        
        % x position of the wave vector
        x = f;
        
        % splitted bands requires additional z offset computations
        isSplitted = (nd == 3) && (floor(iPhase/2) == 1);
        
        % calculate m1 z offset if we are working with 3-D SIM
        if isSplitted
            % beam angle
            %   p = lambda / (2*n*sin(theta)), p is period
            %   sin(theta) = lambda / (2*n*p)
            %              = lambda / (2*n) * f
            %
            % for m1 terms, interference frequency is doubled, hence
            %   sin(theta) = lambda / n * f
            %
            % kz = (1-cos(theta)) * n/lambda
            %    = (1-sqrt(1-sin(theta)^2)) * n/lambda
            %    = n/lambda - sqrt((n/lambda)^2 - f^2)
            y = n/lambda - sqrt((n/lambda)^2 - f^2);
        else
            y = 0;
        end
        
        %% create the impulse functions
        % odd elements are negative terms (offset by m0, so they are even)
        if mod(iPhase, 2) == 0
            x = -x;
        end
        
        % m1 terms of 3-D SIM have splitted bands
        if isSplitted
            % m1 term has splitted bands
            k = [x, x; -y, y];
        else
            k = [x; y];
        end
        
        % turns into the actual wave number
        %   k = 2*pi / p
        %     = 2*pi * f
        % Note: While p is period and f is frequency, we are working 
        % backward from the reciprocal space, hence f is used instead of p.
        k = 2*pi * k;
        
        % splitted bands require additional merging
        % Note: Implicit dimension expansion is used in the wave function.
        if isSplitted
            k = reshape(k, 2, 1, []);
            Im = 1 + real(exp(1i * (vx.*k(2, :, :) + vy.*k(1, :, :))));
            Im = sum(Im, 3);
        else
            Im = 1 + real(exp(1i * (vx.*k(2) + vy.*k(1))));
        end
        
        %DEBUG
        figure(hPre);
        subplot(1, 2, 1);
        imagesc(Im);
            axis image;
            title('R-Space');
        subplot(1, 2, 2);
        imagesc(abs(fftshift(fft2(ifftshift(Im)))).^0.7);
            axis image;
            title('K-Space');
            
        %% apply to system OTF
        TF(:, :, iPhase, iOri) = psf2otf(PSF .* Im, imSz);
        
        %DEBUG
        figure(hPost);
        subplot(1, 2, 1);
        imagesc(PSF);
            axis image;
            title('System PSF');
        subplot(1, 2, 2);
        imagesc(abs(ifftshift(TF(:, :, iPhase, iOri))).^0.7);
            axis image;
            title('Modulated');
    end
end

end
