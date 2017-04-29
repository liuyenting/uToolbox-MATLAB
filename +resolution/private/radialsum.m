function s = radialsum(I, pres)
%RADIALSUM Calculate the radial sum of an image with interpolation.
%
%   S = RADIALSUM(I, ARES) calculate radial sum at different spatial 
%   frequency using angular resolution of ARES. ARES is default to 1 degree
%   if not assigned.

% default pixel resolution to 1px
if nargin == 1
    pres = 1;
end

% find the minimal dimension
[nrows, ncols] = size(I);
if nrows ~= ncols
    warning('resolution:radialsum', ...
            'Not a square matrix, cropped to the minimal square.');
end
% TODO: crop the image

r = 0:floor(nrows/2);
s = zeros([L, 1]);
for i = 1:length(L)
    smpl = radialsmplr(I, r(i), pres);
    s(i) = sum(smpl);
end

end
