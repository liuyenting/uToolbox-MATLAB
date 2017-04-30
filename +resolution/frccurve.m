function [fcfrq, fcraw, fcavg, fcstd] = frccurve(coords, npx, nt, blk)
%FRCCURVE Calculate Fourier ring correlation curve.
%
%   NPX     Super-resolved image size in pixels.
%   NT      N trials to perform the averaging.
%   BLK     N blocks to randomize the dataset.

if npx(1) ~= npx(2)
    error('resolution:frccurve', ...
          'Image size should be a square.');
end

if nargin == 3
    blk = 500;
end

% estimate proper pixel dimensions that can contain all the data
pxsz = estpxsize(coords, npx);

% radial sample size, assuming the dimension are matched
nrs = floor(npx(1)/2)+1;

% start parallel pool
nthread = 24;
if isempty(gcp('nocreate'))
    if nthread > nt
        % shrink the thread pool if not required
        nthread = nt;
    end
    parpool('local', nthread);
end

% ensembeld result
fcraw = zeros([nt, nrs]);
% start the iterations
fprintf('%d tasks, running %d at a time...\n', nt, nthread);
parfor i = 1:nt
    % shuffle the input
    scoords = shuffle(coords, 2, blk);
    
    % bin the data 
    I1 = resolution.binlocal(scoords{1}, npx, pxsz);
    I2 = resolution.binlocal(scoords{2}, npx, pxsz);
    
    % generate the FRC curve
    fcraw(i, :) = loesssmooth(resolution.frc(I1, I2));
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
