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

% normalize the intensity across phases
I = normopint(I, volSz, parms);

for iz = 1:nz
    for iOri = 1:nOri
        for iPhase = 1:nPhase
            disp(num2str(iz));
            % extract volume
            T = I(iz, iOri, iPhase, :, :);
            T = squeeze(T);

            % pad the surrounding sides
            psz = [10, 10];
            T = padarray(T, psz, 0, 'both');
            % RL deconvolution
            T = deconvlucy(T, parms.PSF, parms.PreDeconv);
            % crop the result
            T = T(psz(1)+1:end-psz(1), psz(2)+1:end-psz(2));

            % save the processed volume
            I(iz, iOri, iPhase, :, :) = T;
        end
    end
end

% find kp if it does not exist
if ~isfield(parms, 'Kp') || isempty(parms.Kp)
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
            P = I(:, iOri, iPhase, :, :);
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

    %% find out kp value using the first stack
    
end

%% reconstruction

end
