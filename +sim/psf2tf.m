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

pxSz = parms.PixelSize;

n = parms.RefractiveIndex;
lambda = parms.Wavelength;

%% pre-calculate
OTF = psf2otf(PSF, imSz);

%% pre-allocate
TF = zeros([imSz, nPhase, nOri], 'single');

% first phases require no special treatments
TF(:, :, 1, :) = repmat(OTF, [1, 1, 1, nOri]);

% modulated illumination plane
Im = zeros(imSz, 'single');

% the grid for computation
[vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
% offset to center
midpt = floor(imSz/2)+1;
vx = vx-midpt(1);
vy = vy-midpt(2);
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
            % kz = (1-cos(theta)) * n/lambda
            %    = (1-sqrt(1-sin(theta)^2)) * n/lambda
            %    = n/lambda - sqrt((n/lambda)^2 - (f/2)^2)
            y = n/lambda - sqrt((n/lambda)^2 - (f/2)^2);
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
        
        % Since we have two possible k formation, shift them to ensure we
        % will not collide with the image dimension (X and Y).
        k = reshape(k, 2, 1, []);
        
        % turns into cycle (wave number)
        k = (2*pi) / k;
        
        % Note: Implicit dimension expansion is used here.
        %Im = 1 + cos(vx.*k(1, :, :)) .* cos(vy.*k(2, :, :));
        Im = 1 + cos(vx.*k(1, :, 1));
        
        % merge the result if we are working with splitted bands
        if isSplitted
            Im = sum(Im, 3);
        end
    end
end

end
