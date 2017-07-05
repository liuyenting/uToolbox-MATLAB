function I = normopint(I, sz, parms)
%NORMOPINT Normalize intensities across different illumination setups.
%   Detailed explanation goes here

% extract parameters
nOri = parms.Orientations;
nPhases = parms.Phases;
nz = sz(3);

% find out normalization ratio
ratio = zeros([nz, nOri, nPhases], 'single');
for io = 1:nOri
    for ip = 1:nPhases
        % extract the volume
        S = I(:, io, ip, :, :);
        S = squeeze(S);
        
        % sum along X, Y dimension
        ratio(:, io, ip) = sum(sum(S, 3), 2);
    end
end
% mean along the orientation and phases
mratio = mean(mean(ratio, 3), 2);
% divide all the values with the ratio
ratio = ratio ./ repmat(mratio, [1, nOri, nPhases]);

% normalize
for io = 1:nOri
    for ip = 1:nPhases
        % extract the volume
        N = I(:, io, ip, :, :);
        N = squeeze(N);
        
        % divide the ratio
        r = ratio(:, io, ip);
        N = N ./ repmat(r, [1, sz(1), sz(2)]);
        
        % save to output
        I(:, io, ip, :, :) = N;
    end
end

end

