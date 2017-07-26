function kp = findkp(I, M, parms)
%FINDKP Find the precise values of the pattern wave vector.
%   
%   TBA

%% parameters
volSz = size(I);

imSz = volSz(1:2);

nOri = parms.Orientations;
nPhase = parms.Phases;

%% pre-calculate
% midpoint of current image size
midpt = floor(imSz/2)+1;

hPre = figure('Name', 'XCorr for Kp Estimation', 'NumberTitle', 'off');

%% pre-allocate
% initialize the kp
kp = zeros([2, nPhase-1, nOri], 'single');

%% process
for iOri = 1:nOri
    % convert to frequency space
    D = fftshift(fft2(ifftshift(I(:, :, :, iOri))));
    
    %% retrieve domains
    % flatten the array
    D = reshape(D, [prod(imSz), nPhase]);
    % solve the matrix
    D = (M \ D')';
    % reshape back to original image size
    D = reshape(D, [imSz, nPhase]);
    
    %% find peaks (frequencies)
    for iPhase = 2:nPhase
        %% estimate kp
        X = fxcorr2(D(:, :, 1), D(:, :, iPhase));
        
        % preview
        figure(hPre);
        subplot(nOri, nPhase-1, (iOri-1)*nOri+iPhase-1);
        imagesc(abs(X));
            axis image;
            colormap(gray);
            title(sprintf('O%d P%d', iOri, iPhase-1));
        
        % find the position of the peak
        [~, ind] = max(X(:));
        [y, x] = ind2sub(imSz, ind);
               
        %% parabolic interpolation   
        % position offset from initial guess
        xo = parapeak([X(x-1, y), X(x, y), X(x+1, y)]);
        yo = parapeak([X(x, y-1), X(x, y), X(x, y+1)]);

        % exact position, remove the position offset
        x = x - xo;
        y = y - yo;
        
        % preview
        hold on;
        plot(x, y, 'oy');
        drawnow;
        
        %% from position to shift
        % convert to offset with the center
        x = x - midpt(1);
        y = y - midpt(2);

        kp(:, iPhase-1, iOri) = [x, y];
    end
end

end

function C = fxcorr2(A, B)
%FXCORR2 Fast 2-D cross-correlation.
%
%   C = FXCORR2(A, B) performs cross-correlation upon image A and B. Size
%   of C is the maximum size of A and B on X and Y dimension.
%
%   See also: FFT2, IFFT2, FFTSHIFT, IFFTSHIFT

% % real data only
% if ~isreal(A) || ~isreal(B)
%     error(generatemsgid('InvalidInType'), 'Only real data are allowed.');
% end

% find the region that can cover both A and B
%   size of an image is [nrows (y), ncols (x)]
sz = max(size(A), size(B));

% Since cross-correlation is essentially a convolution, while convolution 
% can be implemented as element-wise multiplication in the reciprocal 
% space, we simply pad the input images A, B to enough size and perform an
% FFT/IFFT, viola!
f1 = fftshift(fft2(ifftshift(A), sz(1), sz(2)));
f2 = fftshift(fft2(ifftshift(B), sz(1), sz(2)));
fx = f1 .* f2;
% C = fftshift(ifft2(ifftshift(fx), 'symmetric'));
C = fftshift(ifft2(ifftshift(fx)));

end

function C = parapeak(vars)
%PARAPEAK Find the vertex of a parabola.
%
%   TBA

%% validate
if length(vars) ~= 3
    error(generatemsgid('InvalidInput'), ...
          'Position variable should be exactly three.');
end

if (vars(1) == vars(2)) && (vars(2) == vars(3))
    warning(generatemsgid('MalformedInput'), ...
            'Unable to solve the input, default to center.');
    
    % directly return the result
    C = 0;
    return;
end

%% compute
% Assuming
%   f(x) = Ax^2 + Bx + C
% and
%   f(-1) = a, f(0) = b, f(1) = c
% so
%   A - B + C = a
%   C = b
%   A + B + C = c
%
% Solving the above lienar equations yields
%   A = (a+c-2*b)/2, B = (c-a)/2
% 
% By deriving f(x), we know the peak is located at
%   x = -B/(2*A)
%     = -(1/2)*(c-a)/(a+c-2*b)
C = -(vars(3)-vars(1))/(vars(1)+vars(3)-2*vars(2)) / 2;

end