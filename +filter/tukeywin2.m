function W = tukeywin2(sz, ratio)
%TUKEYWINN 2-D Tukey window.
%
%   W = TUKEYWIN1(SZ, RATIO) returns a 2-D Tukey (cosine-tapered)
%   window. The RATIO parameter specifies the ratio of the length of taper
%   section to the total length of the window. RATIO should be in the range
%   of [0, 1].
%
%   See also: TUKEYWIN1

Wx = filter.tukeywin1(sz(1), ratio);
Wx = Wx.';

Wy = filter.tukeywin1(sz(2), ratio);

W = Wy .* Wx;

end

