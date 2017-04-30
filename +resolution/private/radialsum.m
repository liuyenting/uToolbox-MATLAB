function s = radialsum(I)
%RADIALSUM Calculate the radial sum of an image with interpolation.
%
%   S = RADIALSUM(I, ARES) calculate radial sum at different spatial 
%   frequency using angular resolution of ARES. ARES is default to 1 degree
%   if not assigned.

% find the minimal dimension
[nrows, ncols] = size(I);
if nrows ~= ncols
    warning('resolution:radialsum', ...
            'Not a square matrix, cropped to the minimal square.');
    I = sqcrop(I);
end

% radius sample location
r = 0:floor(nrows/2);
nr = length(r);
rm = distmap([nrows, ncols]);

% start sampling
s = zeros([nr, 1]);
gI = gpuArray(I);
gR = gpuArray(rm);
for i = 1:nr
    gR = (gR >= r(i)) & (gR < r(i)+1);
    gS = sum(sum(gI .* gR));
    s(i) = gather(gS);
end

end

function r = distmap(sz)
%DISTMAP Generate distance map from the center of the image.

% find the midpoint
midpt = sz / 2;

% generate comparison grid
[xq, yq] = meshgrid(1:sz(1), 1:sz(2));
xq = xq - midpt(1);
yq = yq - midpt(2);

r = sqrt(xq.^2 + yq.^2);

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