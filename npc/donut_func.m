close all;
clearvars;

[vx, vy] = meshgrid(1:32);
v
parm = [211.3974   18.0911   16.8014    3.8418    1.5494];
V = donut(vx, vy, parm);

imagesc(V);
axis image;

function F = donut(vx, vy, parm)

A = parm(1);
xc = parm(2);
yc = parm(3);
R = parm(4);
sigma = parm(5);

[~, r] = cart2pol(vx-xc, vy-yc);

F = A * exp(-(r-R).^2/(2*sigma^2));

end

function F = gauss(vx, parm)

A = parm(1);
c = parm(2);
sigma = parm(5);

F = A * exp(-(vx-c)^2/(2*sigma^2));

end
