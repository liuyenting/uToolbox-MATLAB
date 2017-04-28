function ind = ndsub2ind(sz, sub)
%NDSUB2IND Convert N dimensional subscripts to linear indices.

% cumulative products for the dimension
k = cumprod(sz);
k = [1, k(1:end-1)];

ind = sub * k';

end
