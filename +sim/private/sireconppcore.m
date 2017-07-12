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
    A = filter.tukeywin2(imSz, parms.ApodizeRatio);
    A = single(A);
    
    % positivity constraints
    A(A < 0) = 0;
    
    figure('Name', 'Apodization Function', 'NumberTitle', 'off');
    imagesc(A);
        axis image;
end

% interpolated size
rSz = parms.RetrievalInterpRatio*imSz;

% buffer space for results from the frequency domain, each for the original
% image and the padded image
F = zeros([imSz, nPhase], 'single');
Fp = zeros([rSz, nPhase], 'single');
% interpolated result
Fopt = zeros([rSz, nPhase, nOri], 'single');

% buffer space for the relative matrix (single orientation only)
pr = zeros([rSz, nPhase-1], 'single');
% grids for the relative phase shift matrix
[vx, vy] = meshgrid(1:rSz(1), 1:rSz(2));
vx = vx - (rSz(1)+1);
vy = vy - (rSz(2)+1);

for iOri = 1:nOri
    fprintf('.. o = %d\n', iOri);
    
    %% create relative phase shift matrix
    % revert from position to shift
    shift = kp(:, :, iOri) - imSz.'/2;
    % the ratio in current upsampled dimension
    shift = bsxfun(@rdivide, shift, rSz.');
    % calculate shift in unit spatial frequency, the result is negated
    % since we are trying to shift it back to where it should be
    shift = (2*pi) * (-shift);
    
    for iPhase = 1:nPhase-1
        % fill-in the distance matrix with phase shifts (in unit spatial
        % frequency)
        D = vx*shift(1, iPhase) + vy*shift(2, iPhase);
        % convert to imaginary part in order to apply shift in R-space
        pr(:, :, iPhase) = exp(1i * D);
    end
    
    for iPhase = 1:nPhase
        % extract volume
        T = I(:, :, iPhase, iOri);
        
        % pad the surrounding sides
        T = padarray(T, [psz, psz], 0, 'both');
        % RL deconvolution
        T = deconvlucy(T, parms.PSF, parms.PreDeconv);
        % crop the result
        T = T(psz+1:end-psz, psz+1:end-psz);
        
        % perform Fourier transform without the padded region
        F(:, :, iPhase) = fftshift(fft2(ifftshift(T)));
    end
    
    %% retrieve domains
    % flatten the array
    F = reshape(F, [prod(imSz), iPhase]);
    % solve the matrix
    F = (M \ F')';
    % reshape back to original image size
    F = reshape(F, [imSz, iPhase]);
    
    %% find phases
    % multiply apodization function
    F = F .* A;
    % upsampling to perform FT interpolation in real space
    li = floor((rSz-imSz)/2)+1;
    ui = li+imSz-1;
    Fp(li(1):ui(1), li(2):ui(2), :) = F;
    
    % reference, m_0
    Fref = Fp(:, :, 1);
    Fopt(:, :, 1, iOri) = Fref;
    
    %% search the optimal inital phase
    % unit spatial frequency
    lim = (2*pi) ./ rSz;
    
    options = optimoptions( ...
        'fmincon', ...
        'FiniteDifferenceStepSize', max(lim), ...
        'StepTolerance', 1e-2, ...
        'Display', 'iter-detailed' ...
    );
    p0 = fmincon( ...
        @(x) costfunc(Fref, Fp(:, :, 2:end), x, pr), ...
        [0, 0], ...
        [], [], [], [], ...
        [-pi, -pi], [pi, pi], ...
        [], ...
        options ...
    );

    % round p0 to nearest multiple of unit spatial frequency  
    p0 = round(p0./lim) .* lim;
    fprintf('.... m1=%f, m2=%f\n', p0(1), p0(2));
    
    % save the optimal shifted result
    [~, Fres] = costfunc( ...
        Fref, ...               % m_0
        Fp(:, :, 2:end), ...    % m_i
        p0, ...                 % estimated p0
        pr ...                  % relative phase shifts
    );
    Fopt(:, :, 2:end, iOri) = Fres;
end

%% the actual reconstruction
% sum all the orientations and phases
Fopt = reshape(Fopt, [rSz, nOri*nPhase]);
J = sum(Fopt, 3);
J = fftshift(ifft2(ifftshift(J), 'symmetric'));

%% preview the result
% % show the reconstructed result
% figure('Name', 'Reconstructed', 'NumberTitle', 'off');
% imagesc(J);
%     axis image;
% drawnow;

end

function [S, varargout] = costfunc(Fref, Fp, p0, pr)
%COSTFUNC Cost function to minimize for the phase retrieval algorithm.
%
%   sz: image size
%   Rp: frequency plane
%   p0: initial phase shift
%    p: relative phsae shift, determined by kp

% persistent h;
% 
% if isempty(h) || ~isvalid(h)
%     h = figure('Name', 'Phase Retrieval', 'NumberTitle', 'off');
% end

profile resume;

% interleave the phases since we now have m_i^- and m_i^+
p0 = exp(1i * [-p0; p0]);
% flatten the array for the linear duplication later
p0 = p0(:);
np = length(p0);

Fp = bsxfun(@times, Fp, reshape(p0, [1, 1, np]));

% back to time domain
for ip = 1:np
    Fp(:, :, ip) = fftshift(ifft2(ifftshift(Fp(:, :, ip)))); 
end

% add relative phase shift deduced from kp values (imaginary number in the
% time domain)
Fp = Fp .* pr;

Fp = fftshift(fft2(ifftshift(Fp)));
if nargout == 2
    varargout{1} = Fp;
end

% sum the result to evaluate performance
S = Fref + sum(Fp, 3);

% figure(h);
% imagesc(S);
%     axis image;
% drawnow;

% maximize the function, use negative sign to use fmin* optimizer
S = abs(Fref .* S);
S = -sum(S(:));

% output is required to be double instead of single
S = double(S);

profile off;

end
