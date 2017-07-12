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
    
%     figure('Name', 'Apodization Function', 'NumberTitle', 'off');
%     imagesc(A);
%         axis image;
        
    % reuse for each phase
    A = repmat(A, [1, 1, nPhase]);
    % swap dimension
    A = permute(A, [3, 1, 2]);
end

% interpolated size
rSz = parms.RetrievalInterpRatio * imSz;

% buffer space for results from the frequency domain, each for the original
% image and the padded image
F = zeros([nPhase, imSz], 'single');
Fp = zeros([nPhase, rSz], 'single');
% interpolated result
R = zeros([nOri, nPhase, rSz], 'single');

% buffer space for the relative matrix (single orientation only)
pr = zeros([nPhase-1, rSz], 'single');
% grids for the relative phase shift matrix
[vx, vy] = meshgrid(1:rSz(1), 1:rSz(2));
vx = vx - (rSz(1)+1);
vy = vy - (rSz(2)+1);

for iOri = 1:nOri
    %% create relative phase shift matrix
    % revert from position to shift
    shift = (squeeze(kp(iOri, :, :)) - repmat(imSz/2, [nPhase-1, 1]));
    % the ratio in current upsampled dimension
    shift = shift ./ repmat(rSz, [nPhase-1, 1]);
    % calculate shift in unit spatial frequency, the result is negated
    % since we are trying to shift it back to where it should be
    shift = (2*pi) * (-shift);
    
    for iPhase = 1:nPhase-1
        % fill-in the distance matrix with phase shifts (in unit spatial
        % frequency)
        D = vx*shift(iPhase, 1) + vy*shift(iPhase, 2);
        % convert to imaginary part in order to apply shift in R-space
        pr(iPhase, :, :) = exp(1i * D);
    end
    
    for iPhase = 1:nPhase
        % extract volume
        T = I(iOri, iPhase, :, :);
        T = squeeze(T);
        
        % pad the surrounding sides
        T = padarray(T, [psz, psz], 0, 'both');
        % RL deconvolution
        T = deconvlucy(T, parms.PSF, parms.PreDeconv);
        
        % perform Fourier transform without the padded region
        F(iPhase, :, :) = ...
            fftshift(fft2(ifftshift(T(psz+1:end-psz, psz+1:end-psz))));
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
    % upsampling to perform FT interpolation in real space
    li = floor((rSz-imSz)/2)+1;
    ui = li+imSz-1;
    Fp(:, li(1):ui(1), li(2):ui(2)) = F;
    
    % reference, m_0
    Rref = fftshift(ifft2(ifftshift(squeeze(Fp(1, :, :)))));
    R(iOri, 1, :, :) = abs(Rref);
    
    %% search the optimal inital phase
%     % apply nonlinear optimization
%     problem = struct;
%     problem.objective = @(x) costfunc(Rref, R(2:end, :, :), rSz, x, pr);
%     problem.x0 = [0, 0];
%     problem.solver = 'fminsearch';
%     % additional options
%     options = optimset('Display', 'final', 'TolX', 1e-6);
%     problem.options = options;
%     p0 = fminsearch(problem);

    options = optimoptions( ...
        'fmincon', ...
        'FiniteDifferenceStepSize', pi/10, ...
        'StepTolerance', 1e-2, ...
        'Display', 'iter-detailed' ...
    );
    p0 = fmincon( ...
        @(x) costfunc(Rref, Fp(2:end, :, :), x, pr), ...
        [0, 0], ...
        [], [], [], [], ...
        [-pi, -pi], [pi, pi], ...
        [], ...
        options ...
    );

    fprintf('o=%d, m1=%.2f, m2=%.2f, ratio=%.2f\n', iOri, p0(1), p0(2), p0(2)/p0(1));
    
    % save the optimal shifted result
    [~, Ropt] = costfunc( ...
        Rref, ...               % m_0
        Fp(2:end, :, :), ...    % m_i
        p0, ...                 % estimated p0
        pr ...                  % relative phase shifts
    );
    R(iOri, 2:end, :, :) = Ropt;
end

%% the actual reconstruction
% reshape to ignore differences between orientation and phase dimension
R = reshape(R, [nOri*nPhase, rSz]);
% sum all the orientations and phases
J = squeeze(sum(R, 1));

J = abs(J);

%% preview the result
% show the reconstructed result
figure('Name', 'Reconstructed', 'NumberTitle', 'off');
imagesc(J);
    axis image;
drawnow;

end

function [S, varargout] = costfunc(R0, Rp, p0, pr)
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

np = length(p0);

% interleave the phases and double the element count since we now have 
% m_i^- and p_i^+
p0 = exp(1i * [-p0, p0]);
np = 2*np;
% flatten the array for the linear duplication later
p0 = reshape(p0', [1, np]);

Rp = bsxfun(@times, Rp, reshape(p0, [np, 1, 1]));

for ip = 1:np
    % back to real space
    Rp(ip, :, :) = fftshift(ifft2(ifftshift(squeeze(Rp(ip, :, :))))); 
end

% figure('Name', 'Before', 'NumberTitle', 'off'); 
% subplot(1, 2, 1);
%     imagesc(log(abs(squeeze(Rp(1, :, :))))); 
%         axis image;
%         title('m_1^-');
% subplot(1, 2, 2);
%     imagesc(log(abs(squeeze(Rp(2, :, :)))));
%         axis image;
%         title('m_1^+');
% drawnow;

% add relative phase shift deduced from kp values
Rp = Rp .* pr;

% figure('Name', 'After', 'NumberTitle', 'off'); 
% subplot(1, 2, 1);
%     imagesc(log(abs(squeeze(Rp(1, :, :))))); 
%         axis image;
%         title('m_1^-');
% subplot(1, 2, 2);
%     imagesc(log(abs(squeeze(Rp(2, :, :)))));
%         axis image;
%         title('m_1^+');
% drawnow;

if nargout == 2
    varargout{1} = Rp;
end

% sum the result to evaluate performance
S = R0 + squeeze(sum(Rp, 1));
% ensure we are working with the magnitude instead of imaginary numbers
S = abs(S);

% figure(h);
% imagesc(S);
%     axis image;
% drawnow;

% maximize the function, use negative sign to use fmin* optimizer
% remember to squeeze R0 since it is extracted from a multi-dimension array
S = R0 .* S;
S = -sum(S(:)); % / (sum(R0(:)) * sum(C(:)));

% output is required to be double instead of single
S = double(S);

profile off;

end
