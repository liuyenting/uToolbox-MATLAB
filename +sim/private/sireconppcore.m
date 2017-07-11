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

% buffer space for results from the frequency domain
F = zeros([nPhase, imSz], 'single');
% buffer space for results from the interpolated real space
rSz = parms.RetrievalInterpRatio * imSz;
R = zeros([nPhase, rSz], 'single');
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
    % upsampling to perform FT interpolation in real space
    offset = floor((rSz-imSz)/2)+1;
    R(:, offset(1):offset(1)+imSz(1)-1, offset(2):offset(2)+imSz(2)-1) = F;
    
    % reference frame
    T = fftshift(ifft2(ifftshift(squeeze(R(1, :, :)))));
    R(1, :, :) = abs(T);
    
    % revert from position to shift
    shift = (squeeze(kp(iOri, :, :)) - repmat(imSz/2, [nPhase-1, 1]));
    % the ratio in current upsampled dimension
    shift = shift ./ repmat(rSz, [nPhase-1, 1]);
    % calculate shift in unit spatial frequency, the result is negated
    % since we are trying to shift it back to where it should be
    shift = -(2*pi) * shift;
  
    % create the relative phase shift matrix
    [vx, vy] = meshgrid(1:rSz(1), 1:rSz(2));
    vx = vx - (rSz(1)+1);
    vy = vy - (rSz(2)+1);
    pr = zeros([nPhase-1, rSz], 'single');
    for iPhase = 1:nPhase-1
        % distance matrix
        D = vx*shift(iPhase, 1) + vy*shift(iPhase, 2);
        % convert to imaginary part in order to apply shift in R-space
        pr(iPhase, :, :) = exp(1i * D);
    end
    
    %% search the optimal inital phase
    % search density
    dp = 20;
    % search grid in [0, 2*pi]
    phase = linspace(0, 2*pi, dp); 
    % result
    C = zeros(size(phase));
    
    % types of m_i
    nm = (nPhase-1)/2;
    
    % template for the inital phase
    p0 = zeros([nm, 1], 'single');
    
% %     hc = figure('Name', 'Cost Funtion', 'NumberTitle', 'off');
%     
%     %% search for m1
%     for i = 1:dp
%         %i = floor(dp/4);
%         C(i) = costfunc( ...
%             R(1, :, :), ...     % m_0
%             R(2:3, :, :), ...   % m_i
%             rSz, ...            % image size
%             phase(i), ...          % estimated p0
%             pr(1:2, :, :) ...              % rest of the phases
%         );
%     end
%     [~, ind] = max(C);
%     p0(1) = phase(ind);
%     
% %     figure(hc);
% %     subplot(2, 1, 1);
% %         plot(C);
% %         title('m_1');
% %     drawnow;
%     
%     %% search for m2
%     for i = 1:dp
%         %i = floor(dp/4);
%         C(i) = costfunc( ...
%             R(1, :, :), ...     % m_0
%             R(4:5, :, :), ...   % m_i
%             rSz, ...            % image size
%             phase(i), ...          % estimated p0
%             pr(3:4, :, :) ...              % rest of the phases
%         );
%     end
%     [~, ind] = max(C);
%     p0(2) = phase(ind);
%     
% %     figure(hc);
% %     subplot(2, 1, 2);
% %         plot(C);
% %         title('m_2');
% %     drawnow;
%     
%     disp(p0);
    
    % apply nonlinear optimization
    problem = struct;
    problem.objective = @(x) costfunc(R(1, :, :), R(2:end, :, :), rSz, x, pr);
    problem.x0 = [0, 0];
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = zeros(size(problem.x0));
    problem.ub = (2*pi) * ones(size(problem.x0));
    problem.nonlcon = [];
    problem.solver = 'fmincon';
    % additional options
    options = optimoptions( ...
        'fmincon', ...
        'Display', 'iter-detailed' ...
    );
    problem.options = options;
    p0 = fmincon(problem);

    %% the actual reconstruction
    [~, S] = costfunc( ...
        R(1, :, :), ...     % m_0
        R(2:end, :, :), ...   % m_i
        rSz, ...            % image size
        p0, ...          % estimated p0
        pr ...              % rest of the phases
    );
    %TODO combine the data from all the orientation
    J = S;
    
    figure('Name', 'Reconstructed', 'NumberTitle', 'off');
    imagesc(S);
        axis image;
    drawnow;
end

end

function [S, varargout] = costfunc(R0, Rp, sz, p0, pr)
%COSTFUNC Cost function to minimize for the phase retrieval algorithm.
%
%   sz: image size
%   Rp: frequency plane
%   p0: initial phase shift
%    p: relative phsae shift, determined by kp

persistent h;

if isempty(h) || ~isvalid(h)
    h = figure('Name', 'Phase Retrieval', 'NumberTitle', 'off');
end

np = length(p0);

% interleave the phases since we have m_i^- and p_i^+
p0 = 1i * [-p0, p0];
p0 = exp(p0);
% since we have -/+, element count has to double
np = 2*np;
% flatten the array for the linear duplication later
p0 = reshape(p0', [1, np]);
% duplicate the elements, double the size due to -/+
p0 = repmat(p0, [prod(sz), 1]);
% reshape to the proper matrix size
p0 = reshape(p0, [sz, np]);
% m_i major instead of X/Y major
p0 = permute(p0, [3, 1, 2]);
% %DEBUG preview
% figure('Name', 'Initial Phase Shift', 'NumberTitle', 'off');
% subplot(1, 2, 1);
%     imagesc(mod(imag(squeeze(p0(1, :, :))), 2*pi));
%         axis image;
%         title('m_i^-');
% subplot(1, 2, 2);
%     imagesc(mod(imag(squeeze(p0(2, :, :))), 2*pi));
%         axis image;
%         title('m_i^+');
% drawnow;
% shift the frequency plains to their correct locations
Rp = Rp .* p0;

% revert back to real space
for ip = 1:np
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

% apply modulation pattern
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

% sum the result to evaluate performance
C = sum([R0; Rp], 1);
% ensure we are working with the magnitude instead of imaginary numbers
C = abs(C);

% maximize the function, use negative sign to use fmin* optimizer
% remember to squeeze R0 since it is extracted from a multi-dimension array
S = R0 .* C;
% S = -sum(S(:));
S = sum(S(:));

% % output is required to be double instead of single
% S = double(S);

% squeeze the dimension
C = squeeze(C);

figure(h);
imagesc(C);
    axis image;
drawnow;

if nargout == 2
    varargout{1} = C;
end

end
