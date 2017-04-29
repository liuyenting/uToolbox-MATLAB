function A = radialsmplr(I, r, pres)
%RADIALSMPLR Sampling data in radial direction.
%
%   TBA
%
%   Todo: Adapt to calculate the sampling along an elliptical path.

% ensure the input is square 
[nrows, ncols] = size(I);
if ~ismatrix(I) || (nrows ~= ncols)
    error('resolution:radialsum', 'Not a square matrix.');
end

% find the midpoint
midpt = nrows / 2;

% generate the sampling location
ares = calcares(pres, r);
angles = 0:ares:2*pi;

% generate the Cartesian coordinates
[x, y] = pol2cart(angles, r);

% sample on the matrix
A = interp2(I, x+midpt, y+midpt);

% DEBUG
%printsmpltrace(I, x, y, m);

end

function ares = calcares(pres, r)
%CALCARES Calculate angular resolution from desired pixel resolution.

ares = pres / r;

end

function printsmpltrace(A, x, y, m) %#ok<DEFNU>
%PRINTSMPLTRACE Print the sampled trace on source array.

figure('Name', 'Sampling Trace', 'NumberTitle', 'off'); 

% convert the log-scale magnitude if input array is complex
if ~isreal(A)
    A = 100*log(1+abs(fftshift(A)));
end
imagesc(A);
axis image;

hold on;

% dump the trace
plot(x+midpt, y+midpt);

end
