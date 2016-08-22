clearvars; 
%% Parameters
r = 1;          % um
r_res = 0.1;    % um
px = 64;       % pixels

NA = 1.44;
lambda = 488;   % nm
ni = 1.3;

OVR_SMPL_RATIO = 10;

%% Setup the grids.
% Create the radial grid.
maxRadius = px + 1;
rGrid = 1:(maxRadius*OVR_SMPL_RATIO);

% Resize to fit the over sampling ratio.
rGrid = (rGrid-1) / OVR_SMPL_RATIO;

% Convert to real world unit.
rGrid = rGrid * r_res;

%% Calculate h(x, y, z)
% close all;
figure('Name', 'KD Formula', 'NumberTitle', 'off');
for z = (-4*lambda) : (lambda/8) : (4*lambda)
    % Calculate the Kirchhoff's diffraction formula.
    kdIntegral = KirchhoffDiffraction(NA, lambda, z*1e-3, ni);
    hLut = arrayfun(@kdIntegral.calculate, rGrid);

    % Plot the result.
    plot(rGrid, hLut);
    axis tight;
    xlabel('r (\mum)');
    ylabel('I (a.u.)');
    title(['z = ', num2str(z*1e-3), '(\mum)']);
    axis([-inf inf 0 0.3]);
    drawnow;
end
