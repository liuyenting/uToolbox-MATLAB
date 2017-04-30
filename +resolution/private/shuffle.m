function S = shuffle(A, ns, blk)
%SHUFFLE Shuffle the data to fulfill the FRC requirement.
%
%   SHUFFLE(A, NS) 
%   SHUFFLE(A, NS, BLK) shuffles and splits data array A into NS sets. If 
%   BLK is not defined or invalid, it is default to 1.
%
%   Note:
%   Input array A is assumed to be arranged chronologically.

if nargin == 2
    blk = 1;
end

if blk < 0
    warning('resolution:shuffle', ...
            'Invalid block size, revert to default value 1.');
    blk = 1;
end

ndata = size(A, 1);
% pad the rows to multiple of blocks
rsd = mod(ndata, blk);
if rsd ~= 0
    A = padarray(A, blk-rsd, 0, 'post');
end

% reshape as group of blocks
[nrows, ncols] = size(A);
A = reshape(A, [nrows/blk, ncols*blk]);

% generate the permuted indices
permind = randperm(nrows/blk);
% permute the grouped array
A = A(permind, :);

% reshape back to the original size
A = reshape(A, [nrows, ncols]);

% remove the padded value
A = A(1:ndata, :);

% split the data
S = cell(ns, 1);
for i = 1:ns
    S{i} = A(i:ns:end, :);
end

end

function B = padrows(A, usz, val)
%PADROWS Pad the rows of array to multiple of designated unit size.
%
%   PADROWS(A, USZ) pads rows of matrix A to multiply of unit size USZ.
%   PADROWS(A, USZ, VAL) if VAL is assigned, than it is used as the padded
%   value, otherwise, 0 is used.

if nargin == 2
    val = 0;
end



end
