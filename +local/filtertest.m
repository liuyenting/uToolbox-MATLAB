function filtertest
%FILTERTEST B-spline wavelet filter test program.

close all;
clearvars;

%% Load the stack.
imgStackPath = 'C:\Users\Lattice\Documents\MATLAB\data\local_seg_test\art.tif';
tic;
[rawImgStack, rawStackSize] = tiff.imread(imgStackPath, true);
t = toc;
fprintf(' ** %f seconds to load the stack\n', t);
fprintf(' .. stack has the size of (x, y, z) = (%d, %d, %d)\n', ...
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
tic;

q = 3;
s = 2;

% Level 1.
k1 = convkern(1, q, s);
[imgF1, imgV1] = wavflt(rawImg, k1);
% Show the image.
subplot(1, 5, 2);
imagesc(imgF1); axis tight equal;
sd1 = std(imgF1(:));
title(['F1', ' (SD=', num2str(sd1), ')']);

% Level 2.
k2 = convkern(2, q, s);
imgF2 = wavflt(imgV1, k2);
% Show the image
subplot(1, 5, 3);
imagesc(imgF2); axis tight equal;
sd2 = std(imgF2(:));
title(['F2', ' (SD=', num2str(sd2), ')']);

t = toc;
fprintf(' ** %f seconds to perform wavelet filtering\n', t);

%% Thresholding.
ratio = 1;
thImg = imgF2;
thImg(thImg < ratio*sd1) = 0;
% Show the image
subplot(1, 5, 4);
imagesc(thImg); axis tight equal;
title('Masked');

%% Finding peaks.
tic;

[p, rmImg] = local.FastPeakFind(imgF2, ratio*sd1);
im = imfuse(rawImg, rmImg);
% Show the image
subplot(1, 5, 5);
imagesc(im); axis tight equal;
title('Peaks');

%% Count the peaks.
np = numel(p)/2;
col = p(1:2:end);
row = p(2:2:end);

t = toc;
fprintf(' ** %d peaks are found in %f seconds\n', np, t);

%% Sub-pixel localize.
avh = figure('Name', 'Aerial View', 'NumberTitle', 'off');
set(avh, 'Position', [300, 250, 600, 600]);
zih = figure('Name', 'Zoom In', 'NumberTitle', 'off');
set(zih, 'Position', [1000, 250, 600, 600]);

tic;

boxsz = 11;
boxhsz = floor(boxsz/2);
for i = 1:np
    r = row(i);
    c = col(i);
    
    roi = extractroi(rawImg, [c-boxhsz, c+boxhsz, r-boxhsz, r+boxhsz]);
    if isempty(roi)
        row(i) = -1;
        col(i) = -1;
        continue;
    end
    
%     % Show current location.
%     figure(avh);
%     imagesc(rawImg); axis tight equal manual;
%     hold on;
%     rectangle('Position', [c-boxhsz-0.5, r-boxhsz-0.5, boxsz, boxsz], ...
%               'EdgeColor', 'w');
%     hold off;
% 
%     % Show ROI.
%     figure(zih);
%     imagesc(roi); axis tight equal manual;
%     hold on;
%     plot(1+boxhsz, 1+boxhsz, ...
%          'Marker', 's', 'MarkerEdgeColor', 'r', 'MarkerSize', 20);
    
    % Fit.
    parm = fitgauss(roi, boxsz);
    dr = parm(4);
    dc = parm(2);
    
%     r = (1+boxhsz)+dr;
%     c = (1+boxhsz)+dc;
%     plot(c, r, ...
%          'Marker', 'x', 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
%     hold off;
% 
%     drawnow;
    
    % Write back.
    row(i) = row(i)+dr;
    col(i) = col(i)+dc;
end

t = toc;
fprintf(' ** %f seconds to localize all the local maxima\n', t);

%% Write result to file.
% Convert the unit.
pxRealSz = 102; % [nm]
row = row * pxRealSz;
col = col * pxRealSz;

% Get path string.
[dname, fname, ~] = fileparts(imgStackPath);
csvPath = fullfile(dname, [fname, '.csv']);

% Open the file.
fid = fopen(csvPath, 'w');
fprintf(fid, 'x [nm], y [nm]\n');
% Write all the values.
for i = 1:np
    fprintf(fid, '%.3f, %.3f\n', col(i), row(i));
end
% Close the file.
fclose(fid);

end

function varargout = wavflt(V, k)
%WAVFLT B-spline wavelet filter.

[nrow, ncol] = size(V);

% Pad the array. 
pds = numel(k);
Vp = padarray(V, [pds, pds], 'replicate');

% Convolve column.
for icol = 1:ncol
    Vp(:, icol) = conv(Vp(:, icol), k.', 'same');
end
% Convolve row.
for irow = 1:nrow
    Vp(irow, :) = conv(Vp(irow, :), k, 'same');
end

% Remove the padding.
Vp = Vp(1+pds:end-pds, 1+pds:end-pds);

varargout{1} = V - Vp;
if nargout > 1
    varargout{2} = Vp;
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

function J = extractroi(I, roi)
%EXTRACTROI Extract region of interest from specified image.

cMin = roi(1); 
cMax = roi(2);
rMin = roi(3); 
rMax = roi(4);

% Boundary check.
[ncol, nrow] = size(I);
if (cMin < 1) || (rMin < 1) || (cMax > ncol) || (rMax > nrow)
    J = [];
else
    J = I(rMin:rMax, cMin:cMax);
end

end

function parm = fitgauss(I, boxsz)
%FITGAUSS Fit the Gaussian peak.

persistent problem;

%% Initialize.
if ~exist('problem', 'var') || ...
   ~isstruct(problem) || isempty(fieldnames(problem))
    % Objective function.
    problem.objective = @gaussfunc;
    
    % Initial guess, [A, x0, wx, y0, wy, theta].
    problem.x0 = [-1, 0, boxsz, 0, boxsz];

    % Boundaries.
    boxhsz = floor(boxsz/2);
    problem.lb = [0, -boxhsz, 0, -boxhsz, 0];
    w = boxhsz^2;
    problem.ub = [realmax('double'), boxhsz, w, boxhsz, w];
    
    % Grid.
    problem.xdata = gengrid(boxhsz);
    
    % Solver, default string.
    problem.solver = 'lsqcurvefit';
    % Solver optimization options.
    options = optimoptions('lsqcurvefit', ...
                           'FiniteDifferenceType', 'central', ...
                           'Display', 'off'); 
    problem.options = options;
    
    fprintf(' ** Solver constant parameters initialized!\n');
end

%% Setup parameters for current round.
% Guess amplitude.
problem.x0(1) = max(I(:));
% Raw data.
problem.ydata = double(I);

%% Processing.
[parm, ~, ~, flag] = lsqcurvefit(problem); %#ok<ASGLU>

x = parm(2);
y = parm(4);
    
end

function F = gaussfunc(parm, grid)
%GAUSSFUNC Generate a 2-D Gaussian matrix using specified parameters.
%   F = GAUSSFUNC(PARM, GRID) generates a 2-D Gaussian peak using the PARM
%   configuration on a grid assigned by GRID.
%
%   Note: 
%   PARM is a vector that contains all the possible variables in the form 
%   of
%       PARM = [A, x0, wx, y0, wy, theta]
%   with each element denotes
%       A       Amplitude
%       x0, y0  Center of the peak
%       wx, wy  Spread of the peak

% Extract the parameter.
A = parm(1);
x0 = parm(2);
y0 = parm(4);
wx = parm(3);
wy = parm(5);

% Apply the Gaussian function.
F = A * exp( -( ((grid(:, :, 1)-x0).^2)/(2*wx^2) + ...
                ((grid(:, :, 2)-y0).^2)/(2*wy^2) ) );
            
end

function grid = gengrid(hsz)
%GENGRID Generate fitting grid for the solver.

[x, y] = meshgrid(-hsz:hsz, -hsz:hsz);
grid = zeros(size(x, 1), size(y, 1), 2);
grid(:, :, 1) = x;
grid(:, :, 2) = y;

end