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
        fpath = fullfile(wfDir, fname);
        tiff.imsave(P, fpath);

        % save the result for further processing
        Ip(iOri, iPhase, :, :) = P;
    end
end

%% pseudo widefield image
WF = reshape(Ip, [nOri*nPhase, volSz(1:2)]);
WF = sum(WF, 1);
WF = squeeze(WF);
tiff.imsave(WF, fullfile(wfDir, 'widefield.tif'));

if nargout == 1
    varargout{1} = WF;
end

end

