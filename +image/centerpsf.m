function [B, varargout] = centerpsf(A)
%CENTERPSF Move the PSF to center of the matrix.
%   
%   TBA

% ensure we are working with viable data
nd = ndims(A);
if (nd ~= 2) && (nd ~= 3)
    error('image:centerpsf', ...
          '1-D and N-D (greater than 3) PSF are not supported.');
end

%% remoev noise
T = A;
% normalize to [0, 1]
T = T - min(T(:));
T = T / max(T(:));

average = mean(T(:));
% remove the average value as a simple measure to remove the noise
T = T - average;

% positivity constraints
T(T < 0) = 0;

%% calculate the centroid
% size of an 3-D image is [height, width, depth] -> [y, x, z]
sz = size(T);
w = sum(T(:));
if nd == 2
    [vx, vy] = meshgrid(1:sz(2), 1:sz(1));
    
    % multiply the weight
    cx = vx.*T; 
    cy = vy.*T;
    
    % divide the sum
    cx = sum(cx(:)) / w;
    cy = sum(cy(:)) / w;
    
    % combine the result
    centroid = [cx, cy];
else
    [vx, vy, vz] = meshgrid(1:sz(2), 1:sz(1), 1:sz(3));
    
    % multiply the weight
    cx = vx.*T; 
    cy = vy.*T;
    cz = vz.*T;
    
    % divide the sum
    cx = sum(cx(:)) / sz(2);
    cy = sum(cy(:)) / sz(1);
    cz = sum(cz(:)) / sz(3);
    
    % combine the result
    centroid = [cx, cy, cz];
end

%% shift the PSF
B = A;

% find the offset
origin = (sz-1)/2 + 1;
offset = centroid - origin;
% shift has to be integer
offset = floor(offset);
% circular shift along all the dimensions
B = circshift(B, offset);

end

