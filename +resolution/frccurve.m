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

% target image size
sz = [nd, nd];
% estimate proper pixel dimensions that can contain all the data
pxsz = estpxsize(coords, sz);

%% Frequency
pxPerFreq = radialsum(ones(sz));
% radial sample counts
nrs = length(pxPerFreq);

r = 1:nrs;

if pxsz(1) ~= pxsz(2)
    error('resolution:frccurve', ...
          'Anisotropic scale is not applicable.');
end
% q = r / (nL * pxsz) = r / L
frcFrq = r / (nd*pxsz(1));

%% Correlation
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
        frcSpu(i, :) = spurious(frcFrq, pxPerFreq, frcNum, uncert);
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

function frcSpu = spurious(q, pxCnt, numerator, uncert)

% Note:
% q = pxCnt / L

% % v(q)
% v = numerator ./ pxCnt;
% 
% % H(q)
% H = pdffactor(q, uncert);
% 
% % sinc
% s2 = sinc(pxCnt / 2).^2;
% 
% frcSpu = log(abs(v ./ H ./ s2));
% 
% frcSpu = loess(frcSpu, 10);

% v(q)
v = numerator ./ (2*pi * pxCnt);

% H(q)
H = pdffactor(q, uncert);

% sinc
s2 = sinc(pi * pxCnt).^2;

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

