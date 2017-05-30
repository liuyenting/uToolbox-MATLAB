close all;
clearvars;

I = imread('overview.tif');

% sampling size
smplsz = 32;

%% convert to point list
% image size
sz = size(I);
% total pixels
np = prod(sz);

% total samples
ns = sum(I(:));
% blank list
ptlst = zeros([ns, 2]);

% iterate through list
is = 1;
for ip = 1:np
    % retrieve the value
    sval = I(ip);
    
    % expand
    if sval > 0
        ptlst(is:is+sval-1, :) = ind2sub(sz, ip);
        is = is+sval;
    end
end
assert(ns == is-1);
fprintf('%d points in the list\n', ns);

%% create voronoi diagram
voronoi(ptlst(:, 1), ptlst(:, 2));