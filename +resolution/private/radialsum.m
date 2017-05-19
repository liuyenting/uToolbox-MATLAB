function s = radialsum(I, smplratio, midpt)
%RADIALSUM Calculate the radial sum of an image with interpolation.
%
%   S = RADIALSUM(I, ARES) calculate radial sum at different spatial 
%   frequency using angular resolution of ARES. ARES is default to 1 degree
%   if not assigned.
%   SMPLRATIO implies the sampling ratio, (0, 1], default at 0.5.

% find the minimal dimension
sz = size(I);
if sz(1) ~= sz(2)
    error('resolution:radialsum', ...
            'Not a square matrix, cropped to the minimal square.');
end

if nargin == 1
    % default cutoff at 50% (Nyquist criterion)
    smplratio = 0.5;
else
    if (smplratio <= 0) || (smplratio > 1)
        error('resolution:radialsum', ...
              'Sampling ratio can only range from (0, 1].');
    end
end

%% Find the center of the image.
% radius sample location
r = 1:floor((sz(1)/2) * smplratio);
% find the center
[~, ind] = max(I(:));
%TODO: dirty hack
if ind == 1
    midpt = sz/2;
    xi = midpt(1);
    yi = midpt(2);
else
    [xi, yi] = ind2sub(sz, ind);
end

%% Sample along the circumfrance. 
s = zeros(size(r));
intp = griddedInterpolant(I);
for i = 1:length(s)
    smpl = radialsmplr(intp, sz, r(i), 1, [xi, yi]);
    s(i) = sum(smpl);
end

end
