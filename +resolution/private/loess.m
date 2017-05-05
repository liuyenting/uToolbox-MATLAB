function s = loess(r, nspan)
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
