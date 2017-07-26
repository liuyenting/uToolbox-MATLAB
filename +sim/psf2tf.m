function TF = psf2tf(imSz, PSF, M, parms)
%PSF2TF Convert PSF to transfer functions.
%   
%   TBA

%% parameters
nPhase = parms.Phases;

psfSz = size(PSF);

pxSz = parms.PixelSize;

%% pre-allocate
% domains
D = zeros([imSz, nPhase], 'single');
% the transfer functions
TF = zeros([imSz, nPhase], 'single');

%% retrieve domains (bands separation)
% retrieve the reciprocal space images
for iPhase = 1:nPhase
    D(:, :, iPhase) = fftshift(fft2(ifftshift(PSF(:, :, iPhase)), ...
                                    imSz(2), imSz(1)));
end

% flatten the array
D = reshape(D, [prod(imSz), nPhase]);
% solve the matrix
D = (M \ D')';
% reshape back to original image size
D = reshape(D, [imSz, nPhase]);

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

%% bead mask
%TODO compensate for non-deal point source

%% average and LPF
% cutoff frequency
f = 2*parms.NA / parms.Wavelength;
% cut-off radius in pixel
r = f * (imSz.*pxSz);

midpt = floor(imSz/2)+1;
for iPhase = 1:nPhase
    T = D(:, :, iPhase);
    
    % create radial profile
    % Note: Due to the limitation of radial sampler, only square image is
    % functional for now.
    TF(:, :, iPhase) = statistics.radialmean(T, midpt, floor(r));
end

if parms.Debug
    figure('Name', 'Transfer Functions', 'NumberTitle', 'off');
    for iPhase = 1:nPhase
        % generate title string
        m = floor(iPhase/2);
        if iPhase > 1
            if mod(iPhase, 2) == 0
                s = '^-';
            else
                s = '^+';
            end
        else
            s = '';
        end
        t = sprintf('m_%d%s', m, s);
        
        subplot(1, nPhase, iPhase);
        imagesc(abs(D(:, :, iPhase)).^0.5);
            axis image;
            title(t);
    end
end

end

function err = errfunc(s, Dm, Dp)
%ERRFUNC Cost function for the non-linear fitting of phi.
%
%   ERR = ERRFUNC(S, DM, DP) applies the complex phase shift S (in real) to
%   both minus and plus signs, and generates the differences between them.
%   Output ERR is coerced to the 

err = exp(-1i*s)*Dm - exp(1i*s)*Dp;
% return type is required to be double
err = double(err);

end