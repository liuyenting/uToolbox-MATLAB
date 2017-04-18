close all; clear all;

%% variables
% frames
n = 16;
% grid view size
row = 4;
col = 4;

ccd.XDim = 64;
ccd.YDim = 64;
ccd.PixelSize = 10;

psf.MaxIntensity = 1000;
psf.RadialSigma = 30;
psf.AxialSigma = 250;
psf.SmplNSig = 2;

%% initialize
h_grid = figure('Name', 'Press "Enter" to refresh...', 'NumberTitle', 'off');
h_grid.MenuBar = 'none';
h_grid.ToolBar = 'none';
h_grid.Position = [50, 50, 192*col, 192*row];

nview = row*col;
tmp = zeros([nview, ccd.XDim, ccd.YDim], 'single');

%% rock-n-roll
figure(h_grid);
currkey = 'return';
while strcmp(currkey, 'return')
    for r = 1:row
        for c = 1:col
            v = npcgen([], [30, 60], 150, [20, 40], ccd, 0.2, psf, 2);
            
            subplot(row, col, (r-1)*row + c);
            imagesc(v);
            set(gca, 'YDir', 'normal', 'Visible', 'off');
            axis image;
            colormap(gray);
        end
    end
    
    pause; % wait for a keypress
    currkey = get(gcf, 'CurrentKey'); 
end
close(h_grid);