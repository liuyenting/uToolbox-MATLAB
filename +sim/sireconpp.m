function J = sireconpp(I, parms)
%SIRECONPP Plan-by-plan SI reconstruction.
%
%   TBA

persistent kp;

%% parameters
volSz = size(I);

imSz = volSz(1:2);
%DEBUG process single layer only
if parms.Debug
    offset = 42;
    nz = 1;
else
    nz = volSz(3);
end

%% pre-allocate
% iterate through the layers
J = zeros([parms.RetrievalInterpRatio*imSz, nz], 'single');

%% process
% generate spectral matrix on-the-fly
M = spectramat(parms.Phases, parms.I0, parms.I1);

% probe for the existence of Kp values
if isempty(kp)
    % create projection view along orientations and phases
    Ip = sim.wfproj(I, parms);
    % find the the pattern wave vector for each orientation
    kp = findkp(Ip, M, parms);
end

for iz = 1:nz
    tStart = tic;
    
    fprintf('\tz = %d\n', iz);
    
    % extract the layer
    %DEBUG process single layer only
    if parms.Debug
        L = I(:, :, iz + (offset-1), :, :);
    else
        L = I(:, :, iz, :, :);
    end
    % use reshape instead of squeeze to avoid single orientation get
    % squeezed as well
    sz = size(L);
    sz(3) = [];
    L = reshape(L, sz);
    
    % run the reconstruction on specific layer
    J(:, :, iz) = sireconppcore(L, M, kp, parms);
    
    tElapsed = toc(tStart);
    fprintf('\t%.3fs\n\n', tElapsed);
end

end
