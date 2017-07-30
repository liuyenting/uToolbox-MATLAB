function W = tukeywin2(sz, ratio)
%TUKEYWINN 2-D Tukey window.
%
%   W = TUKEYWIN1(SZ, RATIO) returns a 2-D Tukey (cosine-tapered)
%   window. The RATIO parameter specifies the ratio of the length of taper
%   section to the total length of the window. RATIO should be in the range
%   of [0, 1].
%
%   See also: TUKEYWIN1

[vx, vy] = meshgrid(1:sz(1), 1:sz(2));
vx = vx - sz(1)/2;
vy = vy - sz(2)/2;
% generate the distance matrix
[~, D] = cart2pol(vx, vy);

% radius is limited by the minimal dimension
r = min(sz);

% create the LUT
x = linspace(-r/2, r/2, r);
y = filter.tukeywin1(r, ratio);

% return 0 for out-of-bound values
W = interp1(x, y, D, 'linear', 0);

end

