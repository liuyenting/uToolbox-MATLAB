function [Ip, varargout] = wfproj(I, parms)
%WFPROJ Create widefield projection.
%   
%   TBA

%% parameters
volSz = size(I);

imSz = volSz(1:2);

nOri = parms.Orientations;
nPhase = parms.Phases;

%% create the directory
wfDir = fullfile(parms.DebugPath, 'z_projection');
mkdir(wfDir);

%% create MIP for each orientation and phase
Ip = zeros([imSz, nPhase, nOri], 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % extract the volume
        P = I(:, :, :, iPhase, iOri);

        % MIP along Z axis
        P = max(P, [], 3);

        % save to file
        fName = sprintf('o%d_p%d.tif', iOri, iPhase);
        fPath = fullfile(wfDir, fName);
        tiff.imsave(P, fPath);
        % save the result for further processing
        Ip(:, :, iPhase, iOri) = P;
    end
end

%% pseudo widefield image
WF = reshape(Ip, [imSz, nOri*nPhase]);
WF = sum(WF, 3);
tiff.imsave(WF, fullfile(wfDir, 'widefield.tif'));

if nargout == 1
    varargout{1} = WF;
end

end

