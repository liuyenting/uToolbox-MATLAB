function S = radialsmplr(I, sz, r)
%RADIALSMPLR Sampling data in radial direction.
%
%   TBA
%
%   Todo: Adapt to calculate the sampling along an elliptical path.

% find the midpoint
midpt = sz / 2;

% generate comparison grid
[xq, yq] = meshgrid(1:sz(1), 1:sz(2));
xq = xq - midpt(1);
yq = yq - midpt(2);

rq = sqrt(xq.^2 + yq.^2);
S = I((rq >= r) & (rq < r+1));

% DEBUG
%printsmpltrace(I, xq, yq);

end

function printsmpltrace(A, x, y) %#ok<DEFNU>
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
plot(x, y);

end
