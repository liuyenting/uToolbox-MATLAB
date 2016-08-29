function filtertest
%FILTERTEST B-spline wavelet filter test program.

close all;
clearvars;

%% Load the stack.
imgStackPath = 'C:\Users\Lattice\Documents\MATLAB\data\local_seg_test\dense.tif';
tic;
[rawImgStack, rawStackSize] = tiff.imread(imgStackPath, true);
t = toc;
fprintf(' ** %f seconds to load the stack\n', t);
fprintf(' ** Stack has the size of (x, y, z) = (%d, %d, %d)\n', ...
        rawStackSize(1), rawStackSize(2), rawStackSize(3));

%% Initiate the figure.
h = figure('Name', 'Result', 'NumberTitle', 'off');
set(h, 'Position', [100, 300, 2400, 600]);

%% Select specific image and region.
layer = 42;
rawImg = rawImgStack(:, :, layer);

% Convert to DOUBLE for calculation.
rawImg = double(rawImg);

% Show the image.
subplot(1, 5, 1);
imagesc(rawImg); axis tight equal;
sd0 = std(rawImg(:));
title(['Raw', ' (SD=', num2str(sd0), ')']);

%% Wavelet filter.
% Level 1.
imgF1 = wavflt(rawImg, 1);
% Show the image.
subplot(1, 5, 2);
imagesc(imgF1); axis tight equal;
sd1 = std(imgF1(:));
title(['F1', ' (SD=', num2str(sd1), ')']);

% Level 2.
imgF2 = wavflt(rawImg, 2);
% Show the image
subplot(1, 5, 3);
imagesc(imgF2); axis tight equal;
sd2 = std(imgF2(:));
title(['F2', ' (SD=', num2str(sd2), ')']);

%% Thresholding.
ratio = 1;
thImg = imgF2;
thImg(thImg < ratio*sd1) = 0;
% Show the image
subplot(1, 5, 4);
imagesc(thImg); axis tight equal;
title('Masked');

end

function F = wavflt(I, i)
%WAVFLT B-spline wavelet filter.

q = 3;
s = 2;

if (i <= 0) || (i > 2)
    error('Level has to be 1 or 2.');
elseif i == 1
    k1 = convkern(1, q, s);
    
    V1 = I;
    [nrow, ncol] = size(I);
    % Convolve column.
    for icol = 1:ncol
        V1(:, icol) = conv(V1(:, icol), k1.', 'same');
    end
    % Convolve row.
    for irow = 1:nrow
        V1(irow, :) = conv(V1(irow, :), k1, 'same');
    end
    F = I - V1;
else
    k1 = convkern(1, q, s);
    
    V1 = I;
    [nrow, ncol] = size(I);
    % Convolve column.
    for icol = 1:ncol
        V1(:, icol) = conv(V1(:, icol), k1.', 'same');
    end
    % Convolve row.
    for irow = 1:nrow
        V1(irow, :) = conv(V1(irow, :), k1, 'same');
    end
    
    k2 = convkern(2, q, s);
    
    V2 = V1;
    % Convolve column.
    for icol = 1:ncol
        V2(:, icol) = conv(V2(:, icol), k2.', 'same');
    end
    % Convolve row.
    for irow = 1:nrow
        V2(irow, :) = conv(V2(irow, :), k2, 'same');
    end
    F = V1 - V2;
end
    
end

function k = convkern(i, q, s)
%CONVKERN Generates the convolutional kernel for a level I wavelet filter
%using a B-spline wavelet of Q order, with the scaling factor S.

if (i <= 0) || (i > 2)
    error('Level has to be 1 or 2.');
elseif i == 1
    k = bsbasis(q, s);
else
    % Using "a tour", insert holes as the order increases.
    k0 = bsbasis(q, s);
    k = [];
    for i = 1:numel(k0)
        k = [k, k0(i), 0]; %#ok<AGROW>
    end
    % Remove dangling zero.
    k = k(1:end-1);
end

end

function k = bsbasis(q, s)
%BSBASIS Generates a B-spline basis function of parameter Q, S.

% Maximum range.
l = 2*ceil(q*s/2)-1;
% Indices for the kernel element.
i = 1:l;
% 1-D grid.
x = i - (l+1)/2;

% The kernel.
t = x/s + q/2;
k = bspline(q, t);

% Scale the kernel to have the sum as one.
k = k/sum(k);

end

function A = bspline(q, t)
%BSPLINE Generates the B-spline coefficient.

if q <= 0
    error('Q has to be greater or equal than 1.');
elseif q == 1
    A = (0 <= t) & (t < 1);
else
    A = (t./(q-1)) .* bspline(q-1, t) + ...
        ((q-t)./(q-1)).*bspline(q-1, t-1);
end

end