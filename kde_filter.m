close all;
clearvars;

I = imread('single_npc.tif');

hf = figure('Name', 'Raw Data', 'NumberTitle', 'off');
imagesc(I);
axis image;
colormap(gray);

%% convert to particle list
% number of pixels
np = numel(I);
% number of potential dots
nd = sum(I(:));

% create the list
list = zeros([1, nd]);

id = 1;
% iterate through the pixels
for ip = 1:np
    % repeatance
    n = I(ip);
    
    if n > 0
        % record the position
        list(id:id+n-1) = ip;
        % increase the index
        id = id+n;
    end
end
% records should match the pixel sum
assert(nd == (id-1));

% convert the assignment to coordinates
[y, x] = ind2sub(size(I), list);

% pack as a single variable
data = [x.', y.'];

%% plot the kde 
sz = size(I);
% [bandwidth, density, gx, gy] = kde2d( ...
%     data, ...
%     max(sz), ...    % sampling grid size (squared)
%     [0, 0], sz ...  % sampling range
% ); 
[~, density, gx, gy] = statistics.kde2(I);

figure('Name', 'KDE', 'NumberTitle', 'off');
% plot density on measured grid, image data are transposed
surf(gy, gx, density, 'LineStyle', 'none');
% set the view
view([0, 90]);
set(gca, 'YDir', 'reverse');
axis image;
colormap hot;
colorbar;

%% identify the center
% % calculate the cutoff threshold
% a = mean(density(:));
% s = std(density(:));
% cutoff = a + 2*s;
% 
% mask = density > cutoff;
% 
% se = strel('disk', 5);
% ed = imerode(density, se);
% 
% figure; imagesc(ed); axis image; colormap(gray);

%% k-mean
[vx, vy] = meshgrid(1:sz(1), 1:sz(2));
X = [vx(:), vy(:), density(:)];
[idx, C] = kmeans(X, 1);
C
figure(hf);
hold on;
plot(C(:,1), C(:,2), 'wo');
hold off;