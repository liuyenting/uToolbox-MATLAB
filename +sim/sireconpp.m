function J = sireconpp(I, parms)
%SIRECONPP Plan-by-plan SI reconstruction.
%
%   TBA

% generate spectral matrix on-the-fly
M = spectramat(parms.Phases, parms.I0, parms.I1);

% normalize the intensity across phases
I = 

% find kp if it does not exist
if ~isfield(parms, 'Kp') || isempty(parms.Kp)
    %% create calibration preview
    % create MIP for debug
    Ip = zeros([nOri, nPhase, sz(1:2)], 'single');
    for iOri = 1:nOri
        for iPhase = 1:nPhase
            % extract the volume
            P = I(iOri, iPhase, :, :, :);
            P = squeeze(P);

            % MIP along Z axis
            P = max(P, [], 1);
            P = squeeze(P);

            % save the result for further processing
            Ip(iOri, iPhase, :, :) = P;
        end
    end

    %% find out kp value using the first stack

end

%% reconstruction

end
