function s = radialsum(I, rstep, ares)
%RADIALSUM Calculate the radial sum of an image with interpolation.
%
%   S = RADIALSUM(I, R)
%   S = RADIALSUM(I, R, RES) calculate radial sum at R spatial frequency
%   using resolution of RES. RES is default to 1 degree if not assigned.

% default angular resolution
if nargin == 2
    ares = 1;
end

[smpl, lbl] = radialsmplr(I, rstep, ares);
s = sum(smpl);

%
% DEBUG
%
figure('Name', '[DEBUG] Radial Sum', 'NumberTitle', 'off');

% get the largest element in the image
maxelem = max(I(:)) + 1;

I = 100*log(1+abs(fftshift(I)));
I(lbl) = maxelem;

% plot the original image
imagesc(I);
colormap(gray);
axis image;

end
