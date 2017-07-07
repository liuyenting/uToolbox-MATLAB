function W = tukeywin1(L, ratio)
%TUKEYWINN 1-D Tukey window.
%
%   W = TUKEYWIN1(L, RATIO) returns a L-point Tukey (cosine-tapered)
%   window. The RATIO parameter specifies the ratio of the length of taper
%   section to the total length of the window. RATIO should be in the range
%   of [0, 1].

i = linspace(0, L, L).';

% tapper boundary
%   1 .. tl tl+1 .. th-1 th .. end
%   T .. T   1   ..  1   T  ..  T
tl = 

W = (sin(4*pi*i/L)).^2;
W((i > L/ratio) & (i < L*(1-1/ratio))) = 1;

end

