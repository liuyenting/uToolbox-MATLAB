close all;
clearvars;

I = imread('single_npc_3.tif');

hf = figure('Name', 'Raw Data', 'NumberTitle', 'off');
imagesc(I);
colormap(gray);
axis image;

xr = xlim;
yr = ylim;

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
        n = 1;
        
        % record the position
        list(id:id+n-1) = ip;
        % increase the index
        id = id+n;
    end
end
% % records should match the pixel sum
% assert(nd == (id-1));
list = list(1:id-1);

% convert the assignment to coordinates
[y, x] = ind2sub(size(I), list);

%% find hole
x = x.';
y = y.';
m = [x, y];
[mesh, cylinders] = findTheHoles(m);

%% plot results
% Triangulation
hold on;

triplot(mesh.ConnectivityList, mesh.Points(:,1), mesh.Points(:,2), 'Color', 'yellow');
axis equal;

% Automatically recognized cylinders
for i = 1:length(cylinders)
    lX = x(cylinders{i});
    lY = y(cylinders{i});
    plot([lX; lX(1)],[lY; lY(1)], 'w');
end

%% fit the circle
% [xf, yf, rf] = circfit(x, y);
% 
% figure
% imagesc(I);
% hold on;
% plot(x, y, 'w.');
% hold on;
% rectangle( ...
%     'Position', [xf-rf, yf-rf, rf*2, rf*2], ...
%     'Curvature', [1, 1], ...
%     'LineStyle', '-', ...
%     'EdgeColor', 'r' ...
% );
% title( ...
%     sprintf( ...
%         'Best fit: R = %0.1f; Ctr = (%0.1f,%0.1f)', ...
%         rf, xf, yf ...
%     ) ...
% );
% plot(xf, yf, 'g.');
% xlim([xf-rf-2, xf+rf+2]);
% ylim([yf-rf-2, yf+rf+2]);
% axis equal
% 
% figure(hf);
% hold on;
% plot(xf, yf, 'go');