close all;

sz = 512;

A = linspace(0, 2^16, sz);
B = linspace(0, 2^16, sz/2);

A = repmat(A, [sz, 1]);

B = flip(B);
B = repmat(B, [sz/2, 1]);

A = uint16(round(A));
B = uint16(round(B));

A(129:384, 129:384) = B;

imagesc(A);
axis image;
colormap(gray);

tiff.imsave(A, 'highcontra.tif', true);