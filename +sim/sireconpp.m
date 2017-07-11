function J = sireconpp(I, volSz, parms)
%SIRECONPP Plan-by-plan SI reconstruction.
%
%   TBA

persistent kp;

% extract frequently used parameters
nz = volSz(3);

% generate spectral matrix on-the-fly
M = spectramat(parms.Phases, parms.I0, parms.I1);
% phase shift values
kp = [];

% probe for the existence of Kp values
if isempty(kp)
    % create projection view along orientations and phases
    Ip = sim.wfproj(I, volSz, parms);
    % find the Kp values for each orientations
    kp = findkp(Ip, volSz(1:2), M, parms, true);
end

% iterate through the layers
J = zeros([nz, parms.RetrievalInterpRatio*volSz(1:2)], 'single');
for iz = 1:nz
%     %DEBUG override z
%     iz = floor(nz/2);
%     nz = 1;
    
    fprintf('z = %d\n', iz);
    
    % extract the layer
    L = I(iz, :, :, :, :);
    % use reshape instead of squeeze to avoid single orientation get
    % squeezed as well
    sz = size(L);
    L = reshape(L, sz(2:end));
    
    % run the reconstruction on specific layer
    J(iz, :, :) = sireconppcore(L, volSz(1:2), M, kp, parms);
end

end
