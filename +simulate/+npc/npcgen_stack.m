close all; clear all;

%% variables
% frames
n = 16;

ccd.XDim = 64;
ccd.YDim = 64;
ccd.PixelSize = 10;

psf.MaxIntensity = 1000;
psf.RadialSigma = 30;
psf.AxialSigma = 250;
psf.SmplNSig = 2;

%% initialize
h_pre = figure('Name', 'Preview', 'NumberTitle', 'off');
h_pre.MenuBar = 'none';
h_pre.ToolBar = 'none';
h_pre.Position = [100, 100, 512, 256];

h_cam = figure('Name', 'Camera', 'NumberTitle', 'off');
h_cam.MenuBar = 'none';
h_cam.ToolBar = 'none';
h_cam.Position = [662, 100, 256, 256];

handle.Preview = h_pre;
handle.Camera = h_cam;

stack = zeros([ccd.XDim, ccd.YDim, n], 'single');

%% rock-n-roll
for i = 1:n
    stack(:, :, i) = npcgen(handle, [30, 60], 150, [20, 40], ccd, 0.2, psf, 2);
    disp(['frame ' num2str(i)]);
end

%% save the output
tiff.imsave(stack, 'npc_sim.tif', 1);