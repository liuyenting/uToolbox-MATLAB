function [fcfrq, varargout] = frccurve(coords, nd, blk)
%FRCCURVE Calculate Fourier ring correlation curve.
%
%   NPX     Super-resolved image size in pixels.
%   NT      N trials to perform the averaging.
%   BLK     N blocks to randomize the dataset.

if nargin == 3
    blk = 500;
end

% estimate proper pixel dimensions that can contain all the data
pxsz = estpxsize(coords, [nd, nd]);

% radial sample size, assuming the dimension are matched
nrs = floor(nd/2)+1;

% shuffle the input
coords = shuffle(coords, 2, blk);

% bin the data
sz = [nd, nd];
I1 = resolution.binlocal(coords{1}, sz);
I2 = resolution.binlocal(coords{2}, sz);

% generate the curve
[

% ensembeld result



fcraw = zeros([nt, nrs]);
fcnum = zeros([nt, nrs]);
% start the iterations
fprintf('%d tasks in queue\n', nt);
parfor i = 1:nt
    % shuffle the input
    scoords = shuffle(coords, 2, blk);
    
    % bin the data 
    I1 = resolution.binlocal(scoords{1}, npx, pxsz);
    I2 = resolution.binlocal(scoords{2}, npx, pxsz);
    
    % generate the FRC curve
    [raw, num] = resolution.frc(I1, I2);
    fcraw(i, :) = loesssmooth(raw);
    fcnum(i, :) = num;
end

% calculate the average and error no matter we have complete the
% calculation or not
fcavg = mean(fcraw);
fcstd = std(fcraw);

% generate the frequency scale
fcfrq = 0:nrs-1;
fcfrq = fcfrq / (nrs*pxsz(1));

end

function s = loesssmooth(r, nspan)
%LOESSSMOOTH Use LOESS to smooth the incoming curve.

if nargin == 1
    nspan = 20;
end

nd = length(r);

% smoothing span
sspan = ceil(nd/nspan);
sspan = sspan + (1-mod(sspan, 2));

s = smooth(r, sspan, 'loess');

end
