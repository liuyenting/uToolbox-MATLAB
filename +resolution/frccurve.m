function [frcFrq, frcCrv, varargout] = frccurve(coords, nd, varargin)
%FRCCURVE Calculate Fourier ring correlation curve.
%
%   NPX     Super-resolved image size in pixels.
%   NT      N trials to perform the averaging.
%   BLK     N blocks to randomize the dataset.

p = inputParser;
addOptional(p, 'BlockSize', 500, @isnumeric);
addOptional(p, 'Iterations', 20, @isnumeric);
parse(p, varargin{:});

blk = p.Results.BlockSize;
n = p.Results.Iterations;

% radial sample size, assuming the dimension are matched
nrs = floor(nd/2)+1;

% target image size
sz = [nd, nd];
% estimate proper pixel dimensions that can contain all the data
pxsz = estpxsize(coords, sz);

% generate temporary storage
frcRaw = zeros([n, nrs]);
frcNum = zeros([n, nrs]);

for i = 1:n
    fprintf('... %d / %d\n', i, n);
    
    % suffle the data
    shufCoords = shuffle(coords, 2, blk);
    
    % bin the data
    I1 = resolution.binlocal(shufCoords{1}, sz, pxsz);
    I2 = resolution.binlocal(shufCoords{2}, sz, pxsz);
    
    % generate the curve
    [tmpFrcRaw, frcNum(i, :)] = resolution.frc(I1, I2);
    frcRaw(i, :) = loess(tmpFrcRaw);
end

% generate the frequency scale
frcFrq = 0:nrs-1;
if pxsz(1) ~= pxsz(2)
    error('resolution:frccurve', ...
          'Anisotropic scale is not applicable in FRC.');
else
    pxsz = pxsz(1);
end
frcFrq = frcFrq / (nrs*pxsz);

% post statistics
if n == 1
    frcCrv = frcRaw;
else
    frcCrv = mean(frcRaw);
end

% spurious correction
if nargout >= 3
    
    varargout{1} = frcSpu;
end

end
