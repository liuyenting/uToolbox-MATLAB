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
    % ensure we are not working with imaginary values
    R(1, :, :) = abs(T);

    % create the illumination pattern
    
    
    % apply nonlinear optimization
    problem = struct;
    problem.objective = @(x) costfunc(R, rSz, x, IL);
    problem.x0 = zeros([1, (nPhase-1)/2]);
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = zeros([1, (nPhase-1)/2]);
    problem.ub = repmat(2*pi, [1, (nPhase-1)/2]);
    problem.nonlcon = [];
    problem.solver = 'fmincon';
    problem.options = [];
    p = fmincon(problem);
end

end

function C = costfunc(Rp, sz, p, IL)
%COSTFUNC Cost function to minimize for the phase retrieval algorithm.
%
%   sz: image size
%   Rp: frequency plane
%   p : phases
%   Ip: illumination pattern

np = length(p);

% generate phase shift
P = repmat(p, [prod(sz), 1]);
P = reshape(P, [sz, np]);
P = permute(P, [3, 1, 2]);
% shift the frequency plains to their correct locations
Rp = Rp .* P;

% revert back to real space
for ip = 2:np
    T = fftshift(ifft2(ifftshift(squeeze(Rp(ip, :, :))))); 
    Rp(ip, :, :) = abs(T);
end

% apply modulation pattern
Rp(2:end, :, :) = Rp(2:end, :, :) .* IL;

% sum the result to evaluate performance
C = sum(Rp, 1);
% maximize the function, use negative sign to use fmin* optimizer
C = squeeze(Rp(1, :, :)) .* C;
C = -sum(C(:));

end
