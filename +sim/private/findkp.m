function kp = findkp(I, imSz, parms)
%FINDKP Summary of this function goes here
%   Detailed explanation goes here

% extract frequent use parameters
nOri = parms.Orientations;
nPhase = parms.Phases;
psz = parms.PadSize;

% generate spectral matrix on-the-fly
M = spectramat(parms.Phases, parms.I0, parms.I1);

% buffer space for results from the frequency domain
F = zeros([nPhase, 2*imSz], 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % extract volume
        T = I(iOri, iPhase, :, :);
        T = squeeze(T);
        
        % pad the surrounding sides
        T = padarray(T, [psz, psz], 0, 'both');
        % RL deconvolution
        T = deconvlucy(T, parms.PSF, parms.PreDeconv);
        % crop the result
        T = T(psz+1:end-psz, psz+1:end-psz);
        
        % Fourier transform
        %   Double the transform size to perform interpolation in the
        %   reciprocal space. This behavior is not observed in the ordinary
        %   reconstruction due to possible losses, hence multiplication as
        %   sine wave modulation in spatial domain is preferred.
        T = fftshift(fft2(T, 2*imSz(1), 2*imSz(2)));
        
        % save the result
        F(iPhase, :, :) = T;
    end
    
    % flatten the array
    F = reshape(F, [iPhase, prod(2*imSz)]);
    % solve the matrix
    F = M \ F;
    % reshape back to original image size
    F = reshape(F, [iPhase, 2*imSz]);
    
    % find the m1-/m1+ terms
    tic;
    X = zeros([2, 2*imSz]);
    % process m1-
    F1 = fftshift(fft2(squeeze(F(1, :, :))));
    F2 = fftshift(fft2(squeeze(F(2, :, :))));
    FX = F1 .* F2;
    X(1, :, :) = ifftshift(ifft2(FX));
    % process m1+
    F1 = fftshift(fft2(squeeze(F(1, :, :))));
    F2 = fftshift(fft2(squeeze(F(3, :, :))));
    FX = F1 .* F2;
    X(2, :, :) = ifftshift(ifft2(FX));
    toc;
    
    % preview the result
    figure;
    subplot(1, 2, 1);
        imagesc(abs(squeeze(X(1, :, :))));
        axis image;
    subplot(1, 2, 2);
        imagesc(abs(squeeze(X(2, :, :))));
        axis image;
end

end

