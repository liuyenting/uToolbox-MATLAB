function [Ip, varargout] = wfproj(I, parms)
% WFPROJ Create widefield projection.
%   
%   IP = WFPROJ(I, PARMS) generates the Z projection along all orientations
%   and phases.
%   [IP, WF] = WFPROJ(I, PARMS) WF is the widefield image, which is
%   essentially summing all the orientations and phases together.

%% parameters
volSz = size(I);

imSz = volSz(1:2);

nOri = parms.Orientations;
nPhase = parms.Phases;

%% create MIP for each orientation and phase
Ip = zeros([imSz, nPhase, nOri], 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % extract the volume
        P = I(:, :, :, iPhase, iOri);

        % MIP along Z axis
        P = max(P, [], 3);
        
        % save the result for further processing
        Ip(:, :, iPhase, iOri) = P;
    end
end

%% pseudo widefield image
if (nargout == 1) || parms.Debug
    WF = reshape(Ip, [imSz, nOri*nPhase]);
    WF = sum(WF, 3);
    if nargout == 2
        varargout{1} = WF;
    end
end

%% debug save
if parms.Debug
    wfDir = fullfile(parms.DebugPath, 'z_projection');
    mkdir(wfDir);
    
    % widefield
    tiff.imsave(WF, fullfile(wfDir, 'widefield.tif'));
    % each orientation and phases
    for iOri = 1:nOri
        for iPhase = 1:nPhase
            fName = sprintf('o%d_p%d.tif', iOri, iPhase);
            fPath = fullfile(wfDir, fName);
            tiff.imsave(Ip(:, :, iPhase, iOri), fPath);
        end
    end
end

end

