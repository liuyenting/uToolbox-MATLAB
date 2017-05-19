function [frcFrq, frcCrv, varargout] = frccurve(coords, res, varargin)
%FRCCURVE Calculate Fourier ring correlation curve.
%
%   [X, Y] = FRCCURVE(COORDS, RES) calculates the Fourier correlation plot
%   using the provided coordinates and assigned super-resolved resolution
%   (in nm).
%   [X, Y, S] = FRCCURVE(..., UNCERT) not only calculates the correlation
%   plot, it also calculates the spurious coorelation using the uncertinaty
%   values provided. Notice that matrix dimension of UNCERT must match the 
%   one provided with the COORDS.
%   [...] = FRCCURVE(..., PARMS) can provides additional control over the
%   randomization process. See Additional Parameters sectino for more
%   details.
%
%   Additional Parameters
%   ---------------------
%   'BlockSize'   Control the batch size during randomization, default to 
%                 500 entries per block.
%   'Iterations'  Iterations to average for the final outcome, default to 
%                 20 cycles.
%
%   Note
%   ----
%   It is recommended that super-resolved resolution should be at least
%   four times of the original resolution, e.g. if the original camera
%   pixels are of size 103nm, the RES argument should be less than 25.75nm.

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

% estimate proper image dimensions that can contain all the data
sz = estimsize(coords, res, 'MaxSquare');
fprintf('pixel size = %.2fnm, image size = %dx%d\n', res, sz(1), sz(2));

%% Frequency
pxPerFreq = radialsum(ones(sz));
% radial sample counts
nrs = length(pxPerFreq);

r = 1:nrs;

% q = r / (nL * pxsz) = r / L
L = sz(1)*res;
frcFrq = r / L;

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
    I1 = resolution.binlocal(shufCoords{1}, sz, res);
    I2 = resolution.binlocal(shufCoords{2}, sz, res);
    
    % generate the curve
    [tmpFrcCrv, frcNum] = resolution.frc(I1, I2);
    frcCrv(i, :) = loess(tmpFrcCrv, 20);
    
    if calcSpurious
        frcSpu(i, :) = spurious(frcFrq, L, frcNum, uncert);
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

function frcSpu = spurious(q, L, numerator, uncert)

% Note:
% q = pxCnt / L

% v(q)
v = numerator ./ (q*L);

% H(q)
H = pdffactor(q, uncert);

% sinc
s2 = sinc(q*L/2).^2;

frcSpu = log(abs(v ./ H ./ s2));
frcSpu = loess(frcSpu, 10);

end

function h = pdffactor(q, uncert)

% calculate the statistic among the uncertainty value
uncAvg = mean(uncert);
uncStd = std(uncert);

% common factors
% fac = 1 + 2 * (2*pi * uncStd * q).^2;
fac = 1 + 2 * (uncStd * q).^2;

% power of the exponential term
% pwr = (2*pi * uncAvg * q).^2;
pwr = (uncAvg * q).^2;

h = exp(-pwr) ./ sqrt(fac);

end

