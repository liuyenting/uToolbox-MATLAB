function J = normopint(I, parms)
%NORMOPINT Normalize intensities across different illumination setups.
%   Detailed explanation goes here

% extract parameters
nOri = parms.Orientations;
nPhases = parms.Phases;

% find out normalization ratio
ratio = zeros([nOri, nPhases], 'single');
for io = 1:nOri
    for ip = 1:nPhases
        % extract the volume
        T = I(io, ip, :, :, :);
        
        ratio(io, ip) = sum(T(:));
    end
end
ratio = ratio / mean(ratio(:));

% normalize
for io = 1:nOri
    for ip = 1:nPhases
        % extract the volume
        T = I(io, ip, :, :, :);
        
        ratio(io, ip) = sum(T(:));
    end
end

end

