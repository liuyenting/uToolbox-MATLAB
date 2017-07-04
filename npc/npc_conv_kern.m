close all;
clearvars;

%% load image
I = imread(fullfile(userpath, 'npc_average', 'overview.tif'));

% sampling size
smplsz = 64;

figure('Name', 'Preview', 'NumberTitle', 'off');
imagesc(I);
set(gca, 'CLim', [0, 10]);
axis image;

%% generate donut filter
[vx, vy] = meshgrid(1:smplsz);
parm = [32^2, smplsz/2, smplsz/2, 10, 4];
D = donut(vx, vy, parm);
% rationalize
D = round(D);
% shift downward to apply penalty
Sold = sum(D(:));
for i = 1:max(D(:))
    T = D - i;
    S = sum(T(:));
    if (S <= 0) && (Sold > 0)
        s = i;
        break;
    end
end
fprintf('shift threshold to %d\n', s);
D = D - s;
% rescale
range = [min(D(:)), max(D(:))];
factor = max(abs(range));
D = D / factor;
fprintf('rescale from [%d, %d] to [%.2f, %.2f]\n', range(1), range(2), min(D(:)), max(D(:)));

figure('Name', 'Kernel', 'NumberTitle', 'off');
imagesc(D);
axis image;
colorbar;

%% convolve
J = conv2(I, D, 'same');

figure('Name', 'Filtered', 'NumberTitle', 'off');
imagesc(J);
set(gca, 'CLim', [-200, 200]);
axis image;
colorbar;

% figure;
% J = (J - min(J(:))) / max(J(:));
% imhist(J);

%% viewer
figure('Name', 'Convolved', 'NumberTitle', 'off');

ax1 = subplot(1, 2, 1);
imagesc(I);
set(ax1, 'CLim', [0, 10]);
axis image;

ax2 = subplot(1, 2, 2);
imagesc(J);
set(ax2, 'CLim', [-200, 200]);
axis image;

linkaxes([ax1, ax2], 'xy');

%% threshold
level = graythresh(J);
mask = imbinarize(J, level);

%% filter and find centroids
% fill hole
filled = imfill(mask, 'holes');

% filter by average area size
s = regionprops(filled, 'Area');
areas = cat(1, s.Area);
% get the largest bin
figure('Name', 'ROI Area Size', 'NumberTitle', 'off');
h = histogram(areas);
xlabel('Area Size (px)');
[~, index] = max(h.Values);
range = h.BinEdges(index:index+1);
% filter
selected = bwareafilt(filled, range);

%% viewer
figure('Name', 'Selected Regions', 'NumberTitle', 'off');

ax1 = subplot(1, 2, 1);
imagesc(I);
set(ax1, 'CLim', [0, 10]);
axis image;

ax2 = subplot(1, 2, 2);
imagesc(selected);
axis image;

linkaxes([ax1, ax2], 'xy');

%% find the centroids
s = regionprops(selected, I, 'WeightedCentroid');
centroid = cat(1, s.WeightedCentroid);

%% remove outlier centroid
imsz = size(I);
roisz = [smplsz, smplsz];

clole = floor(centroid-roisz/2);
cupri = clole+roisz-1;
flag = clole < [1, 1] | cupri >= [imsz(2), imsz(1)];
flag = ~any(flag, 2);

% filter
centroid = centroid(flag, :);
ns = size(centroid, 1);

figure('Name', 'Centroid', 'NumberTitle', 'off');

imagesc(I);
set(ax1, 'CLim', [0, 10]);
axis image;

hold on;
plot(centroid(:, 1), centroid(:, 2), 'yx');

%% crop the images
cropped = zeros([roisz, ns], 'single');
figure('Name', 'ROI', 'NumberTitle', 'off');
for is = 1:ns
    coord = centroid(is, :);
    % shift to lower left
    lole = floor(coord - roisz/2);
    % calcualte the upper right
    upri = lole + roisz - 1;
    R = I(lole(2):upri(2), lole(1):upri(1));
    
    imagesc(R);
    axis image;
    colormap(gray);
    drawnow;
    
    cropped(:, :, is) = R;
end

result = sum(cropped, 3);

figure('Name', 'Sum', 'NumberTitle', 'off');
imagesc(result);
colormap(gray);
axis image;

%% helper functions
function F = donut(vx, vy, parm)

A = parm(1);
xc = parm(2);
yc = parm(3);
R = parm(4);
sigma = parm(5);

[~, r] = cart2pol(vx-xc, vy-yc);

F = A * exp(-(r-R).^2/(2*sigma^2));

end