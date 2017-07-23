function TF = psf2tf(imSz, PSF, kp, parms)
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

%% pre-calculate
OTF = psf2otf(PSF, imSz);

hPre = figure('Name', 'Illumination Pattern', 'NumberTitle', 'off');
hPost = figure('Name', 'Transfer Function', 'NumberTitle', 'off');

%% pre-allocate
TF = zeros([imSz, nPhase, nOri], 'single');

% first phases require no special treatments
TF(:, :, 1, :) = repmat(OTF, [1, 1, 1, nOri]);

% the grid for computation
[vx, vy] = meshgrid(1:psfSz(1), 1:psfSz(2));
% convert to real word scale
vx = vx * pxSz(1);
vy = vy * pxSz(2);

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
