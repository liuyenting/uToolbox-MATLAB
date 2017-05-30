close all;
clearvars;

% file
fileName = 'single_npc_2.tif';
filePath = fullfile(userpath, 'npc_average', fileName);

% roi size
roisz = 32;

%% load the file
tic;

I = imread(filePath);
imsz = size(I);

t = toc;
fprintf('%.2fms to load the image\n', t*1e3);

himg = figure('Name', 'Image', 'NumberTitle', 'off');
p0 = himg.Position;

imagesc(I);
colormap(gray);
axis image;

% ensure the type
I = double(I);

%% generate 3-D data set
% only process pixels > 0
flag = I > 0;

% generate the grid
[vx, vy] = meshgrid(1:imsz(2), 1:imsz(1));

% select the points, save in (X, Y, PX)
L = [vx(flag), vy(flag), I(flag)];

%% delaunay graph
tic; 

DT = delaunayTriangulation(L(:, 1), L(:, 2), L(:, 3));

% calculate free boundary surfaces
[FBtri, FBpoints] = freeBoundary(DT);

t = toc;
fprintf('%.2fms to compute DG\n', t*1e3);

hpts = figure('Name', 'Delaunay Triangulation', 'NumberTitle', 'off');

% double the size
p1 = hpts.Position;
hpts.Position = [p1(1:2), 2*p0(3), p0(4)];

ax1 = subplot(1, 2, 1);
surf(vx, vy, I);
axis equal tight;
grid on;

ax2 = subplot(1, 2, 2);
trisurf( ...
    FBtri, ...
    FBpoints(:,1), FBpoints(:,2), FBpoints(:,3), ...
    'FaceColor', 'cyan', ...
    'FaceAlpha', 0.8 ...
);
axis equal tight;
grid on;

% link the properties
lnkdt = linkprop( ...
    [ax1, ax2], ...
    { ...
        'XLim', 'YLim', 'ZLim', ...
        'CameraUpVector', 'CameraPosition', 'CameraTarget' ...
    } ...
);

%% alpha shape
shp = alphaShape(L, 2);

hpts = figure('Name', 'Alpha Shape', 'NumberTitle', 'off');

% double the size
p1 = hpts.Position;
hpts.Position = [p1(1:2), 2*p0(3), p0(4)];

ax1 = subplot(1, 2, 1);
surf(vx, vy, I);
axis equal tight;
grid on;

ax2 = subplot(1, 2, 2);
plot(shp, 'FaceAlpha', 0.1);
axis equal tight;
grid on;

% link the properties
lnkas = linkprop( ...
    [ax1, ax2], ...
    { ...
        'XLim', 'YLim', 'ZLim', ...
        'CameraUpVector', 'CameraPosition', 'CameraTarget' ...
    } ...
);