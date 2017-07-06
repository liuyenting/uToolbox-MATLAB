function J = sireconpp(I, volSz, parms)
%SIRECONPP Plan-by-plan SI reconstruction.
%
%   TBA

% extract frequent use parameters
nOri = parms.Orientations;
nPhase = parms.Phases;
nz = volSz(3);

% generate spectral matrix on-the-fly
M = spectramat(parms.Phases, parms.I0, parms.I1);

% phase shift values
kp = [];

% iterate through the layers
J = zeros([volSz(3), volSz(1:2)], 'single');
for iz = 1:nz
    % extract the layer
    L = I(iz, :, :, :, :);
    % use reshape instead of squeeze to avoid single orientation get
    % squeezed as well
    sz = size(L);
    L = reshape(L, sz(2:end));
    
    if isempty(kp)
        %% create calibration preview
        % create the directory
        calDir = parms.DebugPath;
        status = mkdir(calDir);
        if ~status
            error('sim:sireconpp', ...
                  'Unable to create the calibration directory.');
        end

        % create MIP for debug
        Ip = zeros([nOri, nPhase, volSz(1:2)], 'single');
        for iOri = 1:nOri
            for iPhase = 1:nPhase
                % extract the volume
                P = L(iOri, iPhase, :, :);
                P = squeeze(P);

                % MIP along Z axis
                P = max(P, [], 1);
                P = squeeze(P);

                % save to file
                fname = sprintf('o%d_p%d.tif', iOri, iPhase);
                fpath = fullfile(calDir, fname);
                tiff.imsave(P, fpath);

                % save the result for further processing
                Ip(iOri, iPhase, :, :) = P;
            end
        end

        % save a pseudo widefield image
        WF = reshape(Ip, [iOri*iPhase, volSz(1:2)]);
        WF = sum(WF, 1);
        WF = squeeze(WF);
        tiff.imsave(WF, fullfile(calDir, 'widefield.tif'));

        %% update the Kp value
        kp = findkp(L, volSz(1:2), M, parms);
    end
    
    % run the reconstruction on specific layer
    J(iz, :, :) = sireconppcore(L, volSz(1:2), M, kp, parms);
end

end
