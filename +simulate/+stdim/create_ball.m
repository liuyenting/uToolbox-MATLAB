close all;
clearvars;

I = zeros([128, 128, 128], 'uint16');

cpos = [20, 40, 80];
r = 16;

[vx, vy, vz] = meshgrid(1:128);
R = (vx-cpos(1)).^2 + (vy-cpos(2)).^2 + (vz-cpos(3)).^2;

I(R <= r^2) = 1;