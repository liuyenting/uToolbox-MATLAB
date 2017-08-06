close all;
clearvars;

%% load the data
imageFilePath = fullfile(userpath, 'npc_average', '21000_23000_concatenated stacks.tif');
I = tiff.imread(imageFilePath);

I = double(I);
roiSz = size(I);
roiSz = roiSz(1:2);

%% accumulate the result
K = sum(I, 3);

roiSz = [128, 128];
K = imresize(K, roiSz);

%% fit to a donut
% generate the grid
[vx, vy] = meshgrid(1:roiSz(1), 1:roiSz(2));

parm = [ ...
    max(K(:)), ...          % amplitude
    roiSz/2, ...            % center (in Cartesian)
    max(roiSz/2)/2, ...     % radius of the donut
    1 ...                   % sigma of the Gaussian distribution
];

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxIterations', 1e4, ...
    'MaxFunctionEvaluations', 1e4, ...
    'Display', 'iter');
[coeff, resnorm] = lsqnonlin(@(x)donut_err(vx, vy, K, x), parm, [], [], options);

fprintf('\n.. amplitude = %.2f\n.. (Xc, Yc) = (%.2f, %.2f)\n.. R = %.2f\n.. sigma = %.2f\n', ...
    coeff(1), coeff(2), coeff(3), coeff(4), coeff(5));

coeff = [97.51, 70.56, 65.70, 18, 7];

% overlay the data
L = donut(vx, vy, coeff);

figure('Name', 'Fit Donut', 'NumberTitle', 'off');

subplot(1, 3, 1);
imagesc(K);
title('Raw');
axis image;

subplot(1, 3, 2);
imagesc(L);
title('Fit');
axis image;

subplot(1, 3, 3);
imagesc(K-L);
title('Difference');
axis image;

%% report
r = coeff(4);
s = coeff(5);

outer = r+s;
inner = r-s;

cx = coeff(2);
cy = coeff(3);

figure('Name', 'Result', 'NumberTitle', 'off');

imagesc(K);
colormap(gray);
rectangle('Position', [cx-r, cy-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor', 'yellow');
rectangle('Position', [cx-outer, cy-outer, 2*outer, 2*outer], 'Curvature', [1, 1], 'EdgeColor', 'cyan');
rectangle('Position', [cx-inner, cy-inner, 2*inner, 2*inner], 'Curvature', [1, 1], 'EdgeColor', 'cyan');
axis image;

fprintf('\n.. outer radius = %.2f\n.. midpoint radius = %.2f\n.. inner radius = %.2f\n', ...
    outer, r, inner);

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

function err = donut_err(vx, vy, V, parm)

F = donut(vx, vy, parm);

% calculate RSE
D = (F-V).^2;
err = sqrt(sum(D(:)));

end

