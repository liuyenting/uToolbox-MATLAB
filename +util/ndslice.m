function B = ndslice(A, dim, ind)
%NDSLICE Extract specific index along designed dimension.
%
%   B = ndslice(A, DIM, IND) slices the N dimensional array A along the DIM
%   dimension, and select the index IND.

% number of dimensions
nd = ndims(A);
% size of the array along each dimension
sz = size(A);

% indices are default to all the members
inds = repmat({}, 1, nd);
for id = 1:nd
    if id == dim
        inds(id) = {ind};
    else
        inds(id) = {1:sz(id)};
    end
end

B = A(inds{:});

end

