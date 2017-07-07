function J = sireconppcore(I, imSz, M, kp, parms)
%SIRECONPPCORE Summary of this function goes here
%   Detailed explanation goes here

persistent A;

% extract frequent use parameters
nOri = parms.Orientations;
nPhase = parms.Phases;
psz = parms.PadSize;

% create the apodization function if not exists
if isempty(A)
    [vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
    
    % shift the center
    vx = vx - imSz(1)/2;
    vy = vy - imSz(2)/2;
    % calculate the distance matrix
    A = sqrt(vx.^2 + vy.^2);
    % create the coefficients
    A = pi/(2*parms.ApodizeRatio) * A;
    % apply the cosine mask
    A = cos(A);
    
    % positivity constraints
    A(A < 0) = 0;
    
    A = single(A);
    
    figure('Name', 'Apodization Function', 'NumberTitle', 'off');
    imagesc(A);
        axis image;
        
    % reuse for each phase
    A = repmat(A, [nPhase, imSz]);
end

% buffer space for results from the frequency domain
F = zeros([nPhase, imSz], 'single');
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
        
        % FT
        T = fftshift(fft2(ifftshift(T)));
        
        % save the result
        F(iPhase, :, :) = T;
    end
    
    %% retrieve domains
    % flatten the array
    F = reshape(F, [iPhase, prod(imSz)]);
    % solve the matrix
    F = M \ F;
    % reshape back to original image size
    F = reshape(F, [iPhase, imSz]);
    
    %% find phases
    % multiply apodization function
    F = F .* A;
    
end

end

