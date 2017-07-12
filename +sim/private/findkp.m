function kp = findkp(I, imSz, M, parms, show)
%FINDKP Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    show = false;
end

% extract frequent use parameters
nOri = parms.Orientations;
nPhase = parms.Phases;
psz = parms.PadSize;

% buffer space for results from the frequency domain, x2 upsampling
fSz = parms.KpUpsamplingRatio*imSz;
F = zeros([fSz, nPhase], 'single');
% initialize the kp
kp = zeros([2, nPhase-1, nOri], 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % extract volume
        T = I(:, :, iPhase, iOri);
        
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
        F(:, :, iPhase) = fftshift(fft2(ifftshift(T), fSz(2), fSz(1)));
    end
    
    % preview the result
    if show
        figure( ...
            'Name', 'Reciprocal Space', ...
            'NumberTitle', 'off' ...
        );
        subplot(1, 3, 1);
            imagesc(log(abs(F(:, :, 1))));
            axis image;
            title('m_0');
        subplot(1, 3, 2);
            imagesc(log(abs(F(:, :, 2))));
            axis image;
            title('m_1^-');
        subplot(1, 3, 3);
            imagesc(log(abs(F(:, :, 3))));
            axis image;
            title('m_1^+');
    end
    
    %% retrieve domains
    % flatten the array
    F = reshape(F, [prod(fSz), nPhase]);
    % solve the matrix
    F = (M \ F')';
    % reshape back to original image size
    F = reshape(F, [fSz, nPhase]);
    
    %% find peaks (frequencies)
    % compute the magnitude for peak search
    F = abs(F(:, :, 1:3));
    
    % find the m1-/m1+ terms
    X = zeros([fSz, 2], 'single');
    X(:, :, 1) = fxcorr2(F(:, :, 1), F(:, :, 2));
    X(:, :, 2) = fxcorr2(F(:, :, 1), F(:, :, 3));
    
    % preview the result
    if show
        figure( ...
            'Name', 'xcorr result of m0 and m1 terms', ...
            'NumberTitle', 'off' ...
        );
        subplot(1, 2, 1);
            imagesc(X(:, :, 1));
            axis image;
            title('m_1^-');
        subplot(1, 2, 2);
            imagesc(X(:, :, 2));
            axis image;
            title('m_1^+');
    end
    
    %% calculate kp values
    % find the position of the peak
    X = reshape(X, [prod(fSz), 2]);
    [~, ind] = max(X);
    [y, x] = ind2sub(fSz, ind.');  
    
    % distance toward the origin (center of the image)
    dist = [x, y] - fSz/2;
    %TODO +1 or -1 or 0 here?
    % revert to the original sampling frequency
    dist = dist / parms.KpUpsamplingRatio;
    
    % convert to positions and save them
    kp(:, 1:2, iOri) = (dist + imSz/2).';
    kp(:, 3:4, iOri) = (2*dist + imSz/2).';
end

%% print the result
colname = cell([1, nOri]);
for iOri = 1:nOri
    colname{iOri} = sprintf('Orientation%d', iOri);
end

rowname = cell([1, nPhase-1]);
for iPhase = 2:2:nPhase
    i = iPhase/2;
    rowname{iPhase-1} = sprintf('m%d-', i);
    rowname{iPhase} = sprintf('m%d+', i);
end

kpstr = cell([nPhase-1, nOri]);
for iOri = 1:nOri
    for iPhase = 1:nPhase-1
        kpstr{iPhase, iOri} = sprintf( ...
            '(%.2f, %.2f)', kp(1, iPhase, iOri), kp(2, iPhase, iOri) ...
        );
    end
end

result = array2table( ...
    kpstr, ...
    'VariableNames', colname, ...
    'RowNames', rowname ...
);
disp(result);

end

function C = fxcorr2(A, B)
%FXCORR2 Fast 2-D cross-correlation.
%
%   C = FXCORR2(A, B) performs cross-correlation upon image A and B. Size
%   of C is the maximum size of A and B on X and Y dimension.
%
%   See also: FFT2, IFFT2, FFTSHIFT, IFFTSHIFT

% real data only
if ~isreal(A) || ~isreal(B)
    error('sim:findkp:fxcorr2', 'Only real data are allowed.');
end

% find the region that can cover both A and B
% size of an image is [nrows (y), ncols (x)]
sz = max(size(A), size(B));

% Since cross-correlation is essentially a convolution, while convolution 
% can be implemented as element-wise multiplication in the reciprocal 
% space, we simply pad the input images A, B to enough size and perform an
% FFT/IFFT, viola!
f1 = fftshift(fft2(ifftshift(A), sz(1), sz(2)));
f2 = fftshift(fft2(ifftshift(B), sz(1), sz(2)));
fx = f1 .* f2;
C = fftshift(ifft2(ifftshift(fx), 'symmetric'));

end