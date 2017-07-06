function kp = findkp(I, imSz, M, parms)
%FINDKP Summary of this function goes here
%   Detailed explanation goes here

% extract frequent use parameters
nOri = parms.Orientations;
nPhase = parms.Phases;
psz = parms.PadSize;

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
end

end

