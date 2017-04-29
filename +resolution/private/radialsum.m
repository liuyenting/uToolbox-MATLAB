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
    I = sqcrop(I);
end

r = 0:floor(nrows/2);
nr = length(r);
s = zeros([nr, 1]);
for i = 1:nr
    smpl = radialsmplr(I, r(i), pres);
    s(i) = sum(smpl);
end

end

function Iout = sqcrop(Iin)
%SQCROP Crop input array to square array.

[nrows, ncols] = size(Iin);
if nrows == ncols
    Iout = Iin;
else
    sz = min(nrows, ncols);
    redsz = max(nrows, ncols) - sz;
    
    % padded size is half of the redundant size
    padsz = floor(redsz/2);
    
    if sz == nrows
        Iout = Iin(:, padsz:padsz+sz-1);
    else
        Iout = Iin(padsz:padsz+sz-1, :);
    end
end

end