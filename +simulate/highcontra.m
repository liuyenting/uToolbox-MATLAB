close all;

sz = 501;

% inner size
ssz = floor(sz/2);

A = linspace(0, 2^16, sz);
B = linspace(0, 2^16, ssz);

A = repmat(A, [sz, 1]);

B = flip(B);
B = repmat(B, [ssz, 1]);

A = uint16(round(A));
B = uint16(round(B));

% anchor
c = floor((sz-ssz)/2);
% range
r = c:c+ssz-1;

% paste the square
A(r, r) = B;

imagesc(A);
axis image;
colormap(gray);

tiff.imsave(A, 'highcontra.tif', true);