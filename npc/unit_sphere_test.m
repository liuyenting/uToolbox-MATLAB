close all;
clearvars;

% load the data
load cell4_zcorrected.mat

% center of mass
xm = mean(x); ym = mean(y); zm = mean(z);
% shift the coordinate system
x = x-xm; y = y-ym; z = z-zm;

% calculate the sphere coordinate system
%   th, azimuth
%   phi, elevation
%   r, radius
[th, phi, r] = cart2sph(x, y, z);

% average distance
%   assuming nuclear membrane contains most of the events
rm = mean(r);

% normalize by the mean distance
r = r / rm;
% warp to [0, 2*pi]/[0, pi]
th = th + pi;
phi = phi + pi/2;

% generate the bins (edges only)
tb = linspace(0, 2*pi, 180);
pb = linspace(0, pi, 90);
% generate the unit sphere
[pv, tv] = meshgrid(pb,tb);
xs = sin(pv).*cos(tv); ys = sin(pv).*sin(tv); zs = cos(pv); 

% discretize the data
%   output in index of tb and pb
td = discretize(th, tb);
pd = discretize(phi, pb);

% analyze the distance along the membrane distance
n = numel(r);
% color
c = zeros(size(zs));
% iterate through the index
for i = 1:n
%     % distance
%     d = abs(r(i)-1);
%     if d > 1
%         d = 1;
%     end
%     % strength (inverse proportion)
%     s = 1-d;
    s = 1;
    
    % save the distance result
    ti = td(i); pi = pd(i);
    c(ti, pi) = c(ti, pi) + intensity(i);
end

c = log(c);

% plot
figure('Name', 'Distribution', 'NumberTitle', 'off');
subplot(1, 2, 1);
h1 = surf(xs, ys, zs, c, 'LineStyle', 'none');
    axis image;
    grid off;
    axis off;
    title('North');
    view([0, 30])
subplot(1, 2, 2);
h2 = surf(xs, ys, zs, c, 'LineStyle', 'none');
    axis image;
    grid off;
    axis off;
    title('South');
    view([0, -30])
    
% record video
v = VideoWriter('view_intensity_sum.avi');
v.FrameRate = 36;
v.open();
% rotate the sphere
for i = 0:360-1
    rotate(h1, [0 0 1], 1);
    rotate(h2, [0 0 1], 1);
    drawnow;
    frame = getframe(gcf);
    v.writeVideo(frame);
end
% save
v.close();