close all; clearvars;

%% set data source
dataDir = 'F:\cell4_zp5um_10ms_a2p5s_r7p5s_RFiSHp2aLFCYC+PmeI_clone_2_labeled';

%% load data
imds = imageDatastore(dataDir, 'ReadFcn', @io.fread);

%% parameters
zOffset = -1;
centroid = [];

% generate the grid bins
azBin = linspace(0, 2*pi, 12);
elBin = linspace(0, pi/2, 6); % hemisphere
% unit sphere
[vp, vt] = meshgrid(elBin, azBin);
xgSph = sin(pi/2-vp).*cos(vt); ygSph = sin(pi/2-vp).*sin(vt); zgSph = cos(pi/2-vp);

% generate the bins
azBin = linspace(0, 2*pi, 540);
elBin = linspace(0, pi/2, 135); % hemisphere
% unit sphere
[vp, vt] = meshgrid(elBin, azBin);
xSph = sin(pi/2-vp).*cos(vt); ySph = sin(pi/2-vp).*sin(vt); zSph = cos(pi/2-vp);

% average window
wndSz = 5;
wndBuf = zeros([size(zSph), wndSz]);

hOffset = figure('Name', 'Coverslip Offset', 'NumberTitle', 'off');
hCentroid = figure('Name', 'Centroid', 'NumberTitle', 'off');
hHist = figure('Name', 'Histogram', 'NumberTitle', 'off');
hRawDist = figure('Name', 'Raw Distance Distribution', 'NumberTitle', 'off');

v = VideoWriter('movie_sphere_p.avi');
v.FrameRate = 5;
v.open();

%% run through the images
iFile = 0;
while hasdata(imds)
    I = read(imds);
    I = double(I);
    iFile = iFile+1;
    
    %DEBUG crop X
    I = I(11:end-10, :, :);
    
    %% find the coverslip
    S = sum(sum(I));
    S = squeeze(S);

%     % 1st-order derivative (unit z-step, zero padded)
%     dSdz = diff(S(:));
%     % 2nd-order derivative
%     d2Sdz2 = diff(dSdz(:));
% 
%     [~, zOffset] = max(d2Sdz2);
%     % two derivatives should yields the actual position at i+1
%     zOffset = zOffset+1;
    [~, zOffset] = max(S);

    figure(hOffset);
    % membrane weighting
    plot(S);
    % determined offset
    yl = ylim;
    hold on; line([zOffset, zOffset], yl, 'Color', 'red'); hold off;
    xlabel('Z position');
    
    %% find the centroid
    imSz = size(I);
    % generate the position grid
    [vy, vx, vz] = meshgrid(1:imSz(2), 1:imSz(1), 1:imSz(3));
    
    % total weights, used the squeezed result
    W = sum(S);
    % apply pixel weightings
    cx = vx.*I; cx = sum(cx(:))/W;
    cy = vy.*I; cy = sum(cy(:))/W;
    cz = vz.*I; cz = sum(cz(:))/W;
    % merge the result
    centroid = [cx, cy, cz];
    
    % create projections
    Ixy = sum(I, 3);
    Iyz = squeeze(sum(I, 2));
    Ixz = squeeze(sum(I, 1));
    
    figure(hCentroid);
    subplot(2, 2, 1);
        imagesc(Ixy);
        axis image;
        hold on; plot(cx, cy, 'r*'); hold off;
    subplot(2, 2, 3);
        imagesc(Ixz.');
        axis image;
        hold on; plot(cx, cz, 'r*'); plot(cx, zOffset, 'ro'); hold off;
    subplot(2, 2, 2);
        imagesc(Iyz);
        axis image;
        hold on; plot(cz, cy, 'r*'); plot(zOffset, cy, 'ro'); hold off;
        
    %% convert coordinate
    % offset the position grid w.r.t. the centroid
    vx = vx-cx; vy = vy-cy; 
    % remove z offset
    vz = vz-cz;
    
    % filter the coordinates
    vx = vx(I>0); vy = vy(I>0); vz = vz(I>0);
    % convert to sphere coordinates
    %   TH azimuth, PHI elevation, R radius
    [az, el, r] = cart2sph(vx, vy, vz);
    
    % warp to positive region 
    %   TH [0, 2*pi], PHI [-pi/2, +pi/2]
    az = az + pi;

    % remove negative elevated samples
    az = az(el>=0); r = r(el>=0); el = el(el>=0);
    n = numel(r);
    nFiltered = W-n;
    fprintf('%d samples are filtered from %d samples\n', nFiltered, W);
    
    figure(hHist);
    histogram(r);
    
%     %% filter by p-value
%     mu = mean(r);
%     sigma = std(r);
%     t = (r-mu)/(sigma/sqrt(n));
%     
%     p = 1-tcdf(t, n-1);
%     az = az(p<0.01); el = el(p<0.01); r = r(p<0.01);
%     % update n
%     n = numel(r);
%     fprintf('%d samples left after p-value filter\n', n);
    
    %% calculate regional average distance    
    % discretize the data, save the index of th and phi
    azDis = discretize(az, azBin);
    elDis = discretize(el, elBin);
    
    % regional average distance (baseline)
    d = zeros(size(zSph));
    for i = 1:n
        % index
        azi = azDis(i); eli = elDis(i);
        
        % extract distance
        val = r(i);
%         val = 1;
        
        % add to the result (averaged)
        d(azi, eli) = d(azi, eli) + val;
%         d(azi, eli) = val;
    end
    
    drawnow;
    
    % circular shift the window buffer
    if iFile < wndSz
        wndBuf(:, :, iFile) = d;
        continue;
    end
    % shift
    wndBuf(:, :, 1:wndSz-1) = wndBuf(:, :, 2:wndSz);
    % add new
    wndBuf(:, :, wndSz) = d;
    
    % calculate the dF/F0
    F0 = mean(wndBuf(:));
    Fr = (d-F0)/F0;
    % remove NaN
    Fr(isnan(Fr)) = 0;
    
    %% plot on sphere
    figure(hRawDist);
    surf(xgSph, ygSph, zgSph, 'EdgeColor', 'w', 'FaceColor', 'none');
        axis image; axis off;
    hold on;
    surf(xSph, ySph, zSph, Fr, 'LineStyle', 'none', 'FaceAlpha', 0.75);
%     surf(xSph, ySph, zSph, d, 'LineStyle', 'none', 'FaceAlpha', 0.75);
        axis image; axis off;
        view([0, 30]);
    hold off;
%     imagesc(Fr.');
%         axis image;
%         set(gca, 'YDir', 'normal'); % flip Y-axis back to (0, 0) on lower-left
%         xlabel('Azimuth'); ylabel('Elevation'); 
    title(sprintf('t = %d', iFile));
    % normalize the colorbar
    c = colorbar;
    set(c, 'Units', 'normalized');
    
    drawnow;
        
    % record to video file
    frame = getframe(gcf);
    v.writeVideo(frame);
    
%     %DEBUG early termination
%     if iFile >= 1
%         break;
%     end
end

% save the video 
v.close();
    
