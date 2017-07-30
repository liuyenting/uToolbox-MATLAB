function I = normopint(I, parms)
%NORMOPINT Normalize intensities across different illumination setups.
%   Detailed explanation goes here

%% parameters
volSz = size(I);

nOri = parms.Orientations;
nPhase = parms.Phases;
nz = volSz(3);

%% process
% find out normalization ratio
ratio = zeros([nPhase, nOri, nz], 'single');
for io = 1:nOri
    for ip = 1:nPhase
        % extract the volume
        S = I(:, :, :, ip, io);
        
        % sum along X, Y dimension
        ratio(ip, io, :) = sum(reshape(S, [volSz(1)*volSz(2), nz]));
    end
end
% mean along the orientation and phases
mratio = mean(reshape(ratio, [nOri*nPhase, nz]), 1);
% divide all the values with the ratio
ratio = bsxfun(@rdivide, ratio, reshape(mratio, [1, 1, nz]));

% normalize
ratio = permute(ratio, [3, 1, 2]);
I = bsxfun(@rdivide, I, reshape(ratio, [1, 1, nz, nPhase, nOri]));

end

