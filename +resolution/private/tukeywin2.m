function mask = tukeywin2(npx, nfrac)
%TUKEYWIN2 Generate a 2-D Tukey window.

if ~ismatrix(npx)
    warning('resolution:tukeywin', 'Only 2-D window is generated.');
    npx = npx(1:2);
end

mask_x = tukeywin1(npx(1), nfrac);
mask_x = mask_x';
mask_x = repmat(mask_x, [1, npx(2)]);

mask_y = tukeywin1(npx(2), nfrac);
mask_y = repmat(mask_y, [npx(1), 1]);

mask = mask_x .* mask_y;

end
