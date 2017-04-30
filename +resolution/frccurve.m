function [frc_raw, frc_avg, frc_std] = frccurve(coords, npx, nt, blk)
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

% radial sample size, assuming the dimension are matched
nrs = floor(npx(1)/2)+1;

% ensembeld result
frc_raw = zeros([nt, nrs]);

% create progress bar
h = waitbar(0, '1', 'Name', 'Calculating FRC...', ...
            'CreateCancelBtn', ...
            'setappdata(gcbf, ''canceling'', 1)');
setappdata(h, 'canceling', 0)

% start the iterations
for i = 1:nt
    % check for cancel button press
    if getappdata(h, 'canceling')
        break
    end
    % report current status
    waitbar(i/nt, h, sprintf('%d/%d', i, nt));
    
    % shuffle the input
    scoords = shuffle(coords, 2, blk);
    
    % bin the data 
    I1 = resolution.binlocal(scoords{1}, npx);
    I2 = resolution.binlocal(scoords{2}, npx);
    
    % generate the FRC curve
    frc_raw(i, :) = loesssmooth(resolution.frc(I1, I2));
end
% delete the waitbar, do not try to close it
delete(h);

% calculate the average and error no matter we have complete the
% calculation or not
frc_avg = mean(frc_raw);
frc_std = std(frc_raw);

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
