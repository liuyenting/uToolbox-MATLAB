close all;
clearvars;

I = imread(fullfile(userpath, 'npc_average', 'overview.tif'));

% sampling size
smplsz = 32;

%% preview
figure('Name', 'Preview', 'NumberTitle', 'off');

imagesc(log(I));
axis image;

%% binarize
mask = I > 0;
% fill hole
filled = imfill(mask, 'holes');

% dilate
se = strel('square', 3);
dilated = imdilate(filled, se);

%% statistics
s = regionprops(dilated, 'Area');
areas = cat(1, s.Area);

figure('Name', 'Historgram', 'NumberTitle', 'off');
histogram(areas, 'BinLimits', [0, 3000]);
xlabel('Area (px)');