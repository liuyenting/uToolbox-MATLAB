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
for i = 1:np
    x0 = x(i); 
    y0 = y(i);
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