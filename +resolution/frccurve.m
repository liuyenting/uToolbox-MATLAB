function [frcFrq, frcCrv, varargout] = frccurve(coords, nd, varargin)
%FRCCURVE Calculate Fourier ring correlation curve.
%
%   NPX     Super-resolved image size in pixels.
%   NT      N trials to perform the averaging.
%   BLK     N blocks to randomize the dataset.

p = inputParser;
addOptional(p, 'UncertaintyXY', [], ...
            @(arr) (size(arr, 1) == size(coords, 1)));
addParameter(p, 'BlockSize', 500, @isnumeric);
addParameter(p, 'Iterations', 20, @isnumeric);
parse(p, varargin{:});

blk = p.Results.BlockSize;
n = p.Results.Iterations;
uncert = p.Results.UncertaintyXY;
% set the flag
calcSpurious = ~isempty(uncert);

% radial sample size, assuming the dimension are matched
nrs = floor(nd/2)+1;

% target image size
sz = [nd, nd];
% estimate proper pixel dimensions that can contain all the data
pxsz = estpxsize(coords, sz);

% generate the frequency scale
frcFrq = 0:nrs-1;
if pxsz(1) ~= pxsz(2)
    error('resolution:frccurve', ...
          'Anisotropic scale is not applicable.');
end
frcFrq = frcFrq.' / (nrs*pxsz(1));

frcCrv = zeros([n, nrs]);
if calcSpurious
    frcSpu = zeros([n, nrs]);
end
for i = 1:n
    fprintf('... %d / %d\n', i, n);
    
    % suffle the data
    shufCoords = shuffle(coords, 2, blk);
    
    % bin the data
    I1 = resolution.binlocal(shufCoords{1}, sz, pxsz);
    I2 = resolution.binlocal(shufCoords{2}, sz, pxsz);
    
    % generate the curve
    [tmpFrcCrv, frcNum] = resolution.frc(I1, I2);
    frcCrv(i, :) = loess(tmpFrcCrv, 20);
    
    if calcSpurious
        frcSpu(i, :) = spurious(nd, frcFrq, frcNum, uncert);
    end
end

% post statistics
if n > 1
    frcCrv = mean(frcCrv);
    if calcSpurious
        frcSpu = mean(frcSpu);
    end
end

if calcSpurious
    varargout{1} = frcSpu;
end

end

function frcSpu = spurious(nd, frcFrq, frcNum, uncert)

% calculate the denominator of v(q)
pxcnt = radialsum(ones([nd, nd]));
    
% v(q)
v = frcNum ./ pxcnt;

% H(q)
H = pdffactor(frcFrq, uncert);

% sinc
s2 = sinc(pi * frcFrq * nd).^2;

frcSpu = log(abs(v ./ H ./ s2));

frcSpu = loess(frcSpu, 10);

end

function h = pdffactor(q, uncert)

% calculate the statistic among the uncertainty value
uncAvg = mean(uncert);
uncStd = std(uncert);

% common factors
fac = 1 + 2 * (2*pi * uncStd * q).^2;

% power of the exponential term
pwr = (2*pi * uncAvg * q).^2;

h = exp(-pwr) ./ sqrt(fac);

end

