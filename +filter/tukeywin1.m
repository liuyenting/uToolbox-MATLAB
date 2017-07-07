function W = tukeywin1(L, ratio)
%TUKEYWINN 1-D Tukey window.
%
%   W = TUKEYWIN1(L, RATIO) returns a L-point Tukey (cosine-tapered)
%   window. The RATIO parameter specifies the ratio of the length of taper
%   section to the total length of the window. RATIO should be in the range
%   of [0, 1].
%
%   Reference
%   ---------
%   https://www.mathworks.com/help/signal/ref/tukeywin.html

% ensure we are working in the valid range
if ~isscalar(L)
    error(generatemsgid('InvalidSize'), 'Invalid dimension size.');
end
if (ratio < 0) || (ratio > 1)
    error(generatemsgid('OutOfBound'), 'Ratio should fall in [0, 1].');
end

% linear grid for the window
i = linspace(0, 1, L).';

% symmetric, divide the ratio in half
period = ratio/2;

% tapper boundary
%   1 .. tl tl+1 .. th-1 th .. end
%   T .. T   1   ..  1   T  ..  T
tl = floor(period*(L-1))+1;
th = L-tl+1;

% create the window in segments
W = [ ...
    ((1 + cos(pi/period * (i(1:tl) - period))) / 2); ...
    ones(th-tl-1, 1); ...
    ((1 + cos(pi/period * (i(th:end) - 1 + period))) / 2)
];

end

