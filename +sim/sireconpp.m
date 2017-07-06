function J = sireconpp(I, volSz, parms)
%SIRECONPP Plan-by-plan SI reconstruction.
%
%   TBA

% extract frequently used parameters
nz = volSz(3);

% generate spectral matrix on-the-fly
M = spectramat(parms.Phases, parms.I0, parms.I1);

% phase shift values
kp = [];

% iterate through the layers
J = zeros([nz, volSz(1:2)], 'single');
for iz = 1:nz
    % extract the layer
    L = I(iz, :, :, :, :);
    % use reshape instead of squeeze to avoid single orientation get
    % squeezed as well
    sz = size(L);
    L = reshape(L, sz(2:end));
    
    % probe for Kp value if we haven't done so
    if isempty(kp)
        kp = findkp(L, volSz(1:2), M, parms);
    end
    
    % run the reconstruction on specific layer
    J(iz, :, :) = sireconppcore(L, volSz(1:2), M, kp, parms);
end

end
