function local_test

close all;
clearvars;

%% Load the stack.
imgStackPath = 'C:\Users\Lattice\Documents\MATLAB\data\local_seg_test\art.tif';
tic;
[rawImgStack, rawStackSize] = tiff.imread(imgStackPath, true);
t = toc;
fprintf(' ** %f seconds to load the stack\n', t);
fprintf(' ** Stack has the size of (x, y, z) = (%d, %d, %d)\n', ...
        rawStackSize(1), rawStackSize(2), rawStackSize(3));
    
%% Select specific image and region.
layer = 42;
rawImg = rawImgStack(:, :, layer);
showimg(rawImg, 'Raw Image');

%% Filter the image.
% Using DWT to denoise the image.
tic;
fltImg = dwtdenoise(rawImg);
t = toc;
fprintf(' ** %f seconds to filter the image\n', t);
showimg(fltImg, 'Filtered Image');

%% Find the peaks.
tic;
p = local.FastPeakFind(fltImg);
t = toc;
np = size(p, 1)/2;
fprintf(' ** %f seconds to find %d peaks in the view\n', t, np);

x = p(1:2:end);
y = p(2:2:end);

%% Open CSV file for saving.
fid = fopen('C:\Users\Lattice\Documents\MATLAB\data\local_seg_test\art_my.csv', 'w');
fprintf(fid, 'x [nm], y [nm], sigma x [nm], sigma y [nm], intensity [au]\n');
%% Fit the result.
avh = figure('Name', 'Aerial View', 'NumberTitle', 'off');
set(avh, 'Position', [500, 100, 600, 800]);
zih = figure('Name', 'Region of Interest', 'NumberTitle', 'off');
set(zih, 'Position', [1200, 200, 600, 600]);

% Size of the bounding box.
boxPxLen = 3;   % [px]
pxLen = 102;    % [nm]
for i = 1:np
    tic;
    
    x0 = x(i); 
    y0 = y(i);
    
    %% Extract ROI.
    % Calculate the bounding box, [Xmin, Xmax, Ymin, Ymax].
    hs = floor(boxPxLen/2);
    bBox = [x0-hs, x0+hs, y0-hs, y0+hs];
    % Extract the ROI.
    roiImg = extractroi(fltImg, bBox);
    % Skip this particle if valid ROI is too small.
    if isempty(roiImg)
        continue;
    end
    
    %% Fitting.
    % Fit the ROI by LSE.
    fitParm = fitgauss(roiImg, x0, y0, boxPxLen);
   
    t = toc;
    fprintf(' ** %f seconds to process\n', t);
    
    %% Plot the result.
    % Aerial view.
    figure(avh);
    imagesc(fltImg);
    title(['Particle ', num2str(i)]);
    hold on;
    % Note: First point is included, therefore the N-0.5f
    rBox = [x0-hs-0.5, y0-hs-0.5, boxPxLen, boxPxLen];
    rectangle('Position', rBox, 'EdgeColor', 'w');
    axis tight equal;
    
    % Zoom in.
    figure(zih);
    imagesc(roiImg);
    hold on;
    plot(x0-rBox(1), y0-rBox(2), 'Marker', 'x', 'MarkerEdgeColor', 'w');
    xf = fitParm(2);
    yf = fitParm(4);
    plot(xf-rBox(1), yf-rBox(2), 'Marker', 's', 'MarkerEdgeColor', 'w');
    axis tight equal;
    hold off;
    xf = xf * 102;
    yf = yf * 102;
    ti = sprintf('(%.3f, %.3f)', xf, yf);
    title(ti);
    
    drawnow;
    
    %% Write to file.
    fprintf(fid, '%.3f, %.3f, %.3f, %.3f, %.3f\n', ...
            xf, yf, fitParm(3)*102, fitParm(5)*102, fitParm(1));
end
fclose(fid);

end

function h = showimg(I, txt)
%SHOWIMG Show the image and resize the figure.

if nargin == 1
    h = figure;
else
    h = figure('Name', txt, 'NumberTitle', 'off');
end
imagesc(I);
axis tight equal;
matchsize(h, I);

end

function matchsize(H, I)
%MATCHSIZE Match the figure to size of its image.

scnSize = get(0, 'ScreenSize');
scnSize = scnSize(3:4);

imgSize = size(I);
figPos = (scnSize - imgSize)/2;
set(H, 'Position', [figPos, imgSize]);

end

function J = dwtdenoise(I)
%DWTDENOISE Denoise the image using Discrete Wavelet Transform.

%% Settings.
% Wavelet family.
wname = 'sym6';
% Decomposition level. 
level = 2;
% Thresholding type.
thrType = 's';
% Threshold limits.
thrSet = [ ...
    14.882939856329342,     14.882939856329342; ...
    61.577863646972794,     57.679169036040989; ...
    14.882939856329342,     14.882939856329342  ...
];

%% Processing.
% Convert to double for processing.
I = double(I);

% Decomposition.
[coefs, sizes] = wavedec2(I, level, wname);
% Denoise.
[J, ~, ~] = wdencmp('lvd', coefs, sizes, wname, level, thrSet, thrType);

%% Output.
% Convert back the result.
J = round(J);
J = uint16(J);

end

function J = extractroi(I, roi)
%EXTRACTROI Extract region of interest from specified image.

xMin = roi(1); 
xMax = roi(2);
yMin = roi(3); 
yMax = roi(4);

% Boundary check.
[ny, nx] = size(I);
if (xMin < 1) || (yMin < 1) || (xMax > nx) || (yMax > ny)
    J = [];
else
    J = I(yMin:yMax, xMin:xMax);
end

end

function parm = fitgauss(I, x0, y0, roiLen)
%FITGAUSS Fit the Gaussian peak.

persistent problem;

%% Initialize.
if ~exist('problem', 'var') || ...
   ~isstruct(problem) || isempty(fieldnames(problem))
    % Objective function.
    problem.objective = @gaussfunc;
    
    % Initial guess, [A, x0, wx, y0, wy, theta].
    problem.x0 = [-1, -1, roiLen, -1, roiLen, 0];

    % Boundaries.
    problem.lb = [0, -1, 0, -1, 0, -pi/4];
    w = (roiLen/2)^2;
    problem.ub = [realmax('double'), -1, w, -1, w, pi/4];
    
    % Solver, default string.
    problem.solver = 'lsqcurvefit';
    % Solver optimization options.
    options = optimoptions('lsqcurvefit', ...
                           'FiniteDifferenceType', 'central'); 
    problem.options = options;
    
    fprintf(' ** Solver constant parameters initialized!\n');
end

%% Setup parameters for current round.
problem.xdata = gengrid(x0, y0, roiLen);
% Guess amplitude.
problem.x0(1) = max(I(:));
% Initial center.
problem.x0(2) = x0;
problem.x0(4) = y0;
% Lower bound.
problem.lb(2) = x0 - roiLen/2;
problem.lb(4) = y0 - roiLen/2;
% Upper bound.
problem.ub(2) = x0 + roiLen/2;
problem.ub(4) = y0 + roiLen/2;
% Raw data.
problem.ydata = double(I);

%% Processing.
[parm, ~, ~, flag] = lsqcurvefit(problem);

% Interpret the flag.
switch(flag)
    case 1
        msg = 'converged';
    case {2, 3, 4}
        msg = 'in tolerance range';
    otherwise
        msg = 'failed';
end
x = parm(2);
y = parm(4);
fprintf(' .. (%d, %d) -> (%.10f, %.10f), %s\n', x0, y0, x, y, msg);
    
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
%       theta   Rotation angle of the peak

% Extract the parameter.
A = parm(1);
x0 = parm(2);
y0 = parm(4);
wx = parm(3);
wy = parm(5);
theta = parm(6);

% Rotate the data points.
rotGrid(:, :, 1) = grid(:, :, 1)*cos(theta) - grid(:, :, 2)*sin(theta);
rotGrid(:, :, 2) = grid(:, :, 1)*sin(theta) + grid(:, :, 2)*cos(theta);

% Rotate the estimated center.
xRot = x0*cos(theta) - y0*sin(theta);
yRot = x0*cos(theta) + y0*sin(theta);

% Apply the Gaussian function.
F = A * exp( -( ((rotGrid(:, :, 1)-xRot).^2)/(2*wx^2) + ...
                ((rotGrid(:, :, 2)-yRot).^2)/(2*wy^2) ) );
            
end

function grid = gengrid(x0, y0, sz)
%GENGRID Generate fitting grid for the solver.

hs = floor(sz/2);
[x, y] = meshgrid(x0-hs:x0+hs, y0-hs:y0+hs);
grid = zeros(size(x, 1), size(y, 1), 2);
grid(:, :, 1) = x;
grid(:, :, 2) = y;

end