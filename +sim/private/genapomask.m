function M = genapomask(imSz, ratio)
%GENAPOMASK Generate the apodization mask.
%   
%   M = GENAPOMASK(IMSZ, RATIO) generates the apodization mask M of size
%   IMSZ. Taper region is controlled by RATIO, which should falls under 
%   [0, 1].

% for now, we default to use Tukey window
M = filter.tukeywin2(imSz, ratio);

end

