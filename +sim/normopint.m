function I = normopint(I, sz, parms)
%NORMOPINT Normalize intensities across different illumination setups.
%   Detailed explanation goes here

% extract parameters
nOri = parms.Orientations;
nPhases = parms.Phases;
nz = sz(3);

% find out normalization ratio
ratio = zeros([nPhases, nOri, nz], 'single');
for io = 1:nOri
    for ip = 1:nPhases
        % extract the volume
        S = I(:, :, :, ip, io);
        
        % sum along X, Y dimension
        ratio(ip, io, :) = sum(reshape(S, [sz(1)*sz(2), sz(3)]));
    end
end
% mean along the orientation and phases
mratio = mean(reshape(ratio, [nOri*nPhases, nz]), 1);
% divide all the values with the ratio
ratio = bsxfun(@rdivide, ratio, reshape(mratio, [1, 1, nz]));

% normalize
ratio = permute(ratio, [3, 1, 2]);
I = bsxfun(@rdivide, I, reshape(ratio, [1, 1, nz, nPhases, nOri]));

end

