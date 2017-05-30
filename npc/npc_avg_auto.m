close all;
clearvars;

I = imread('overview.tif');

% sampling size
smplsz = 32;

%% preview
hfig = figure('Name', 'Preview', 'NumberTitle', 'off');
hfig.Visible = 'off';

imagesc(I);
axis image;
hold on;

%% binarize
mask = I > 0;
% fill hole
filled = imfill(mask, 'holes');

% dilate
se = strel('square', 3);
dilated = imdilate(filled, se);

% keep object in range
selected = bwareafilt(dilated, [250, 300]);

%% retrieve the centroids
s = regionprops(selected, I, 'WeightedCentroid');
centroids = cat(1, s.WeightedCentroid);

%% remove outlier centroid
imsz = size(I);
roisz = [smplsz, smplsz];

flag = centroids(:, [1, 2]) <= roisz | centroids(:, [1, 2]) >= imsz-roisz;
flag = ~any(flag, 2);

% filter
centroids = centroids(flag, :);
ns = size(centroids, 1);

%% preview
msg = sprintf('%d samples selected', ns);
fprintf([msg, '\n']);

plot(centroids(:,1), centroids(:,2), 'yx');
hold off;
title(msg);
hfig.Visible = 'on';

%% crop the images
cropped = zeros([roisz, ns], 'uint16');
for is = 1:ns
    c = centroids(is, :);
    % offset to lower left
    c = floor(c-roisz/2);

    cropped(:, :, is) = uint16(I(c(2):c(2)+smplsz-1, c(1):c(1)+smplsz-1));
end

% save to stack
tiff.imsave(cropped, 'overview_roi.tif', true);

%% align and sum
figure('Name', 'Aligning', 'NumberTitle', 'off');

% coordinate position
[vx, vy] = meshgrid(1:roisz(1), 1:roisz(2));
% template region
result = zeros(roisz);
for is = 1:ns
    J = cropped(:, :, is);
    
    imagesc(J);
    axis image;
    colormap(gray);
    
    % create the list
    lst = im2lst(J);
    
%     [~, d, ~, ~] = statistics.kde2(J);
    [~, d, ~, ~] = kde2d( ...
        lst, ...
        max(roisz), ...    % sampling grid size (squared)
        [0, 0], roisz ...  % sampling range
    ); 

    % compose the result
    X = [vx(:), vy(:), d(:)];
    [~, C] = kmeans(X, 1);
    C = C(1:2)
    
    hold on;
    plot(C(1), C(2), 'yo');
    hold off;
    
    drawnow;
    
    % find shift
    shift = C - roisz/2;
    % align the image
    J = imtranslate(J, -shift, 'OutputView', 'same');
    
    result = result + single(J);
end

figure('Name', 'Sum', 'NumberTitle', 'off');
imagesc(result);
colormap(gray);
axis image;

tiff.imsave(uint16(result), 'overview_sum.tif', true);

function lst = im2lst(im)

% size
sz = size(im);
% number of pixels
np = prod(sz);
% number of potential dots
nd = sum(im(:));

% create the list
lst = zeros([1, nd]);

id = 1;
% iterate through the pixels
for ip = 1:np
    % repeatance
    n = im(ip);
    
    if n > 0
        % record the position
        lst(id:id+n-1) = ip;
        % increase the index
        id = id+n;
    end
end
% records should match the pixel sum
assert(nd == (id-1));

% convert the assignment to coordinates
[y, x] = ind2sub(sz, lst);

% pack as a single variable
lst = [x.', y.'];

end

% %% convert to point list
% % image size
% sz = size(I);
% % total pixels
% np = prod(sz);
% 
% % total samples
% ns = sum(I(:) > 0);
% % blank list
% ptlst = zeros([ns, 2]);
% 
% % iterate through list
% is = 1;
% for ip = 1:np
%     % retrieve the value
%     sval = I(ip);
%     
%     % expand
%     if sval > 0
%         [y, x] = ind2sub(sz, ip);
%         ptlst(is, :) = [x, y];
%         is = is+1;
%     end
% end
% assert(ns == is-1);
% fprintf('%d points in the list\n', ns);
% 
% plot(ptlst(:, 1), ptlst(:, 2), 'r.', 'MarkerSize', 0.1);
% axis equal tight;
% hold on;
% 
% xr = xlim;
% yr = ylim;
% 
% hfig.Visible = 'on';
% 
% %% create voronoi diagram
% data = datasample(ptlst, 200);
% ns = size(data, 1);
% 
% plot(data(:, 1), data(:, 2), 'yx');
% fprintf('probing with %d points\n', ns);
% 
% [vx,vy] = voronoi(data(:, 1), data(:, 2));
% 
% plot(vx, vy, 'b-');
% xlim(xr);
% ylim(yr);
