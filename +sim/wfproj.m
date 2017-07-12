function [Ip, varargout] = wfproj(I, volSz, parms)
%WFPROJ Create widefield projection.
%   Detailed explanation goes here

%% create the directory
wfDir = fullfile(parms.DebugPath, 'widefield');
if exist(wfDir, 'dir') == 7
    content = dir(wfDir);
    if numel(content) > 2
        warning('sim:wfproj', ...
                'Widefield output directory is not empty, overwrite.');
        util.rmcontent(wfDir);
    end
else 
    % create the directory
    status = mkdir(wfDir);
    if ~status
        error('sim:wfproj', ...
              'Unable to create the widefield output directory.');
    end
end

% extract commonly used parameters
nOri = parms.Orientations;
nPhase = parms.Phases;

%% create MIP for each orientation and phase
Ip = zeros([volSz(1:2), nPhase, nOri], 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % extract the volume
        P = I(:, :, :, iPhase, iOri);

        % MIP along Z axis
        P = max(P, [], 3);

        % save to file
        fname = sprintf('o%d_p%d.tif', iOri, iPhase);
        fpath = fullfile(wfDir, fname);
        tiff.imsave(P, fpath);

        % save the result for further processing
        Ip(:, :, iPhase, iOri) = P;
    end
end

%% pseudo widefield image
WF = reshape(Ip, [volSz(1:2), nOri*nPhase]);
WF = sum(WF, 3);
tiff.imsave(WF, fullfile(wfDir, 'widefield.tif'));

if nargout == 1
    varargout{1} = WF;
end

end

