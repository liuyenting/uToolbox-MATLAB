clear all; close all; %#ok<CLALL>

%% loading the data
fprintf('\n -- loading the data --\n');
coords = dlmread(fullfile(userpath, 'subarea3_frc.dat'));

% resolution [dx, dy, dz] in nm
%pxsize = [103, 103, 1000];
pxsize = [103, 103, 1000];
% magnification
mag = 10;

fprintf(' %d samples loaded\n', size(coords, 1));

%% prepare the data set
fprintf('\n -- prepare the data set --\n');

% leave only XY values
coords = coords(:, 1:2);
% offset back to the origin and drop the t-axis
coords = offsetorigin(coords);

% estimate the output size
%[npx, pxsize] = estsize(coords, pxsize, mag);
%fprintf(' x=%d, y=%d, z=%d\n', npx(1), npx(2), npx(3));
%fprintf(' x=%d, y=%d\n', npx(1), npx(2));

%% calculate FRC
fprintf('\n -- calculate FRC --\n');

h = figure('Name', 'FRC resolution', 'NumberTitle', 'off');
hold on;

n = 20;
frc_ens = [];   % ensembled result
for i = 1:n   
    fprintf('%d/%d\n', i, n);
    
    scoords = shuffle(coords, 2, 500);
    
    % bin the data
    I1 = resolution.binlocal(scoords{1}, [1024, 1024]);
    I2 = resolution.binlocal(scoords{2}, [1024, 1024]);
    
    % generate the FRC curve
    frc_cor = resolution.frc(I1, I2);
    
    %nr = radialsum(ones(size(I1)));

    figure(h);
    frc_frq = 0:length(frc_cor)-1;
    frc_frq = frc_frq / size(I1, 1);
    plot(frc_frq, frc_cor);
    %axis([frc_frq(1), frc_frq(end), -0.5, 1]);
    axis tight;
    xlabel('Spatial Frequency (nm^{-1})');
    
    drawnow;
    
    if isempty(frc_ens)
        frc_ens = frc_res;
    else 
        frc_ens = (frc_ens+frc_res) / 2;
    end
end



%frc_frq = 0:length(frc_ens)-1;
%frc_frq = frc_frq / size(I0p, 1);
%plot(frc_frq, frc_ens, 'LineWidth', 2);
%axis([frc_frq(1), frc_frq(end), -0.5, 1]);
%xlabel('Spatial Frequency (nm^{-1})');
hold off;
