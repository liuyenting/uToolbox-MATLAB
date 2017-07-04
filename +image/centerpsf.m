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

% remove the average value as a simple measure to remove the noise
average = mean(A(:));
T = A - average;

% calculate the centroid
sz = size(T);
nelem = prod(sz);
if nd == 2
    [vx, vy] = meshgrid(1:sz(2), 1:sz(1));
    
    % multiply the weight
    cx = vx.*T; 
    cy = vy.*T;
    
    % divide the sum
    cx = sum(cx(:)) / nelem;
    cy = sum(cy(:)) / nelem;
    
    % combine the result
    centroid = [cx, cy];
else
    [vx, vy, vz] = meshgrid(1:sz(2), 1:sz(1), 1:sz(3));
    
    % multiply the weight
    cx = vx.*T; 
    cy = vy.*T;
    cz = vz.*T;
    
    % divide the sum
    cx = sum(cx(:)) / nelem;
    cy = sum(cy(:)) / nelem;
    cz = sum(cz(:)) / nelem;
    
    % combine the result
    centroid = [cx, cy, cz];
end

centroid

B = A;

end

