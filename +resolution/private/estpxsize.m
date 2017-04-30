function pxsz = estpxsize(coords, sz)
%ESTPXSIZE Estimate the pixel size of the binned image.
%
%   PXSZ = ESTPXSIZE(COORDS, SZ) estimates the generated pixel size 
%   according to the maxmium size of the coordinate range.

maxcoord = max(coords);
pxsz = maxcoord ./ (sz-1);
pxsz = repmat(max(pxsz(:)), size(sz));

end
