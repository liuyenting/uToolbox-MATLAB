function local_test

close all;

%% Load the stack.
imgStackPath = 'C:\Users\Lattice\Documents\MATLAB\data\local_seg_test\sparse.tif';
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
      
%% Fit the result.
% Size of the bounding box.
boxSize = 7;

boxHalfSize = floor(boxSize/2);
for i = 1:np
    x0 = x(i); 
    y0 = y(i);
    
    % Calculate the bounding box range, [Xmin, Xmax, Ymin, Ymax].
    bBox = [x0 - boxHalfSize, x0 + boxHalfSize, ...
            y0 - boxHalfSize, y0 + boxHalfSize];
    % Extract the ROI.
    roiImg = extractroi(fltImg, bBox);
    % Skip this particle if valid ROI is too small.
    if isempty(roiImg)
        continue;
    end
    
    % Fit the ROI by LSE.
    [xf, yf] = 
end

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
[nx, ny] = size(I);
if (xMin < 1) || (yMin < 1) || (xMax > nx) || (yMax > ny)
    J = [];
else
    J = I(xMin:xMax, yMin:yMax);
end

end

function parm = fitgauss(x, y)
%FITGAUSS Fit the Gaussian peak.

%% Default parameter.
A_def = 1;
x0_def = 

end

function F = gaussfunc(parm, data)
%GAUSSFUNC Generate a 2-D Gaussian matrix using specified parameters.
%
%   Note: parm = [A, x0, wx, y0, wy, theta]
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
dataRot(:, :, 1) = data(:, :, 1)*cos(theta) - data(:, :, 2)*sin(theta);
dataRot(:, :, 2) = data(:, :, 1)*sin(theta) + data(:, :, 2)*cos(theta);

% Rotate the estimated center.
xRot = x0*cos(theta) - y0*sin(theta);
yRot = x0*cos(theta) + y0*sin(theta);

% Apply the Gaussian function.
F = A * exp( -( ((dataRot(:, :, 1)-xRot).^2)/(2*wx^2) + ...
                ((dataRot(:, :, 2)-yRot).^2)/(2*wy^2) ) );
            
end
