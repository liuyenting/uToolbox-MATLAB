function kp = findkp(I, M, parms)
%FINDKP Find the precise values of the pattern wave vector.
%   
%   TBA

%% parameters
volSz = size(I);

imSz = volSz(1:2);

nOri = parms.Orientations;
nPhase = parms.Phases;

%% pre-allocate
% initialize the kp
kp = zeros([2, nPhase-1, nOri], 'single');

%% process
for iOri = 1:nOri
    % convert to frequency space
    D = fftshift(fft2(ifftshift(I(:, :, :, iOri)), imSz(2), imSz(1)));
    
    %% retrieve domains
    % flatten the array
    D = reshape(D, [prod(imSz), nPhase]);
    % solve the matrix
    D = (M \ D')';
    % reshape back to original image size
    D = reshape(D, [imSz, nPhase]);
    
    %% find peaks (frequencies)
    % compute the magnitude for peak search
    D = abs(D);
    
    % solve each term
    X = zeros([imSz, 2], 'single');
    for iPhase = 2:2:nPhase
        %% standard FFT cross-correlation
        % find the m_i -/+ terms
        X(:, :, 1) = fxcorr2(D(:, :, 1), D(:, :, iPhase));
        X(:, :, 2) = fxcorr2(D(:, :, 1), D(:, :, iPhase+1));
        
        % preview the result
        if parms.Debug
            figure( ...
                'Name', 'xcorr result of m_- and m_+ terms', ...
                'NumberTitle', 'off' ...
            );
            subplot(1, 2, 1);
                imagesc(X(:, :, 1));
                axis image;
                title('m_i^-');
            subplot(1, 2, 2);
                imagesc(X(:, :, 2));
                axis image;
                title('m_i^+');
        end

        % find the position of the peak
        X = reshape(X, [prod(imSz), 2]);
        [~, ind] = max(X);
        [y, x] = ind2sub(imSz, ind);  
        
        %% parabolic interpolation
        X = reshape(X, [imSz, 2]);
        
        % iterate through the -/+ terms
        for i = 1:2
            x0 = x(i);
            y0 = y(i); 
            
            % position offset from initial guess
            xo = parapeak([X(x0-1, y0, i), X(x0, y0, i), X(x0+1, y0, i)]);
            yo = parapeak([X(x0, y0-1, i), X(x0, y0, i), X(x0, y0+1, i)]);
            
            % exact position, remove the position offset
            x(i) = x0 - xo;
            y(i) = y0 - yo;
        end
        
        %% from position to shift
        % distance toward the origin (center of the image)
        midpt = floor(imSz/2)+1;
        dist = [x; y] - midpt.';
        
        % save them
        kp(:, iPhase-1:iPhase, iOri) = dist;
    end
end

%% print the result
% generate column and row labels
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

% generate (x, y) coordinate pair
kpstr = cell([nPhase-1, nOri]);
for iOri = 1:nOri
    for iPhase = 1:nPhase-1
        kpstr{iPhase, iOri} = sprintf( ...
            '(%.2f, %.2f)', kp(1, iPhase, iOri), kp(2, iPhase, iOri) ...
        );
    end
end

% create the table and print-out
result = array2table(kpstr, ...
                     'VariableNames', colname, 'RowNames', rowname);
fprintf('\n');
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
    error(generatemsgid('InvalidInType'), 'Only real data are allowed.');
end

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
C = fftshift(ifft2(ifftshift(fx), 'symmetric'));

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