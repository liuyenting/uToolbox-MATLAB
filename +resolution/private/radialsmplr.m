function A = radialsmplr(ip, sz, r, pres, midpt)
%RADIALSMPLR Sampling data in radial direction.
%
%   TBA
%
%   Todo
%   ----
%   Adapt to calculate the sampling along an elliptical path.


% generate the sampling location
ares = calcares(pres, r);
angles = 0:ares:2*pi;

% generate the Cartesian coordinates
[x, y] = pol2cart(angles, r);

% sample on the matrix
A = ip(x+midpt(1), y+midpt(2));

% DEBUG
%printsmpltrace(I, x, y, midpt);

end

function ares = calcares(pres, r)
%CALCARES Calculate angular resolution from desired pixel resolution.

ares = pres / r;

end

function printsmpltrace(A, x, y, midpt) %#ok<DEFNU>
%PRINTSMPLTRACE Print the sampled trace on source array.

figure('Name', '[DEBUG] Sampling Trace', 'NumberTitle', 'off'); 

% convert the log-scale magnitude if input array is complex
if ~isreal(A)
    A = 100*log(1+abs(fftshift(A)));
end
imagesc(A);
axis image;

hold on;

% dump the trace
plot(x+midpt(1), y+midpt(2));

end