clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'subarea3_frc.dat'));

% resolution [dx, dy, dz] in nm
pxsize = [103, 103, 1000];
% magnification
mag = 10;

fprintf(' %d samples loaded\n', size(coords, 1));

%% prepare the data set
fprintf('\n -- prepare the data set --\n');

% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

% estimate the output size
[npx, pxsize] = estsize(coords, pxsize, mag);
fprintf(' x=%d, y=%d, z=%d\n', npx(1), npx(2), npx(3));

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

h = figure('Name', 'FRC resolution', 'NumberTitle', 'off');
hold on;

n = 20;
frc_ens = [];
for i = 1:n   
    fprintf('%d/%d\n', i, n);
    
    % permuted indices
    permInd = randperm(size(coords, 1));
    % permute the data
    coords = coords(permInd, :);
    
    % bin the data
    I0 = resolution.binlocal(coords(2:2:end, :), npx, pxsize);
    I1 = resolution.binlocal(coords(1:2:end, :), npx, pxsize);
    % generate Z projection
    I0p = sum(I0, 3);
    I1p = sum(I1, 3);
    
    % generate the FRC curve
    frc_res = resolution.frc(I0p, I1p, npx, pxsize);

    figure(h);
    frc_frq = 0:length(frc_res)-1;
    frc_frq = frc_frq / size(I0p, 1);
    plot(frc_frq, frc_res);
    axis([frc_frq(1), frc_frq(end), -0.5, 1]);
    xlabel('Spatial Frequency (nm^{-1})');
    
    drawnow;
    
    if isempty(frc_ens)
        frc_ens = frc_res;
    else 
        frc_ens = (frc_ens+frc_res) / 2;
    end
end

frc_frq = 0:length(frc_ens)-1;
frc_frq = frc_frq / size(I0p, 1);
plot(frc_frq, frc_ens, 'LineWidth', 2);
axis([frc_frq(1), frc_frq(end), -0.5, 1]);
xlabel('Spatial Frequency (nm^{-1})');
hold off;
