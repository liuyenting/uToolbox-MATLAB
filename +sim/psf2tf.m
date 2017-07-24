function TF = psf2tf(imSz, PSF, M, kp, parms)
%PSF2TF Convert PSF to transfer functions.
%   
%   TBA

%% parameters
nPhase = parms.Phases;

psfSz = size(PSF);

pxSz = parms.PixelSize;

lambda = parms.Wavelength;

%% pre-allocate
% domains
D = zeros([imSz, nPhase], 'single');
% the transfer functions
TF = zeros([imSz, nPhase], 'single');

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
if parms.Debug && false
    dispFlag = 'iter-detailed';
else
    dispFlag = 'none';
end
stepSz = max(2*pi ./ imSz);
for iPhase = 2:2:nPhase
    Dm = D(:, :, iPhase);
    Dp = D(:, :, iPhase+1);
   
    options = optimoptions(@lsqnonlin, ...
                           'Display', dispFlag, ...
                           'FiniteDifferenceType', 'central', ...
                           'FiniteDifferenceStepSize', stepSz);
    s = lsqnonlin(@(s)(errfunc(s, Dm, Dp)), pi, [], [], options);
    fprintf('\tk0 (m%d) = %.4f\n', iPhase/2, s);
    
    % apply the cancellations
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
    % functional for now.
    TF(:, :, iPhase) = statistics.radialmean(T, midpt, floor(r));
end

end
