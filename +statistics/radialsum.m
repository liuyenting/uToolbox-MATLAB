function s = radialsum(I, smplratio)
%RADIALSUM Calculate the radial sum of an image with interpolation.
%
%   TBA
%
%   TODO use cutoff pixel position instead of ratio
%
%   See also RADIALSMPLR.

%% parameters
sz = size(I);
if sz(1) ~= sz(2)
    error(generatemsgid('NotSquare'), 'Not a square matrix.');
end

if nargin == 1
    % default cutoff at 50% (Nyquist criterion)
    smplratio = 0.5;
else
    if (smplratio <= 0) || (smplratio > 1)
        error(generatemsgid('InvalidInput'), ...
              'Sampling ratio can only range from (0, 1].');
    end
end

%% find the center of the image.
% radius sample location
r = 1:floor((sz(1)/2) * smplratio);
% find the center
[Imax, ind] = max(I(:));
Imin = min(I(:));
% if min and max are the same, there is no possible a peak exists
if Imin == Imax
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
    smpl = radsmplr(intp, [xi, yi], r(i), 1);
    s(i) = sum(smpl);
end

end
