function I = radialmean(I, midpt, rc)
%RADIALMEAN Calculate the radial mean of an image.
%
%   J = RADIALMEAN(I, MIDPT, RC) generates a radially averaged image J of I
%   using MIDPT as the center and RC as the maximum sampling radius.
%   Out-of-range sample points are treated as 0.
%
%   See also RADIALSMPLR.

%% parameters
imSz = size(I);
if imSz(1) ~= imSz(2)
    error(generatemsgid('NotSquare'), 'Not a square matrix.');
end

if nargin == 2
    % default cutoff at half way of the image (Nyquist criterion)
    rc = floor(imSz/2);
else
    if (length(rc) > 1) 
        if any(rc(rc ~= rc(1)))
            error(generatemsgid('InvalidInput'), ...
                  'Anisotropic sampling radius is not supported yet.');
        else
            rc = rc(1);
        end
    end
end

%% pre-allocate
% radial profile
rp = zeros(size(rc));

% radial size
r = linspace(0, rc-1, rc);

%% process
% out-of-bound values are assigned as NaN
ip = griddedInterpolant(I, 'spline', 'none');
for i = 1:rc
    % always sample with pixel density of 1
    smpl = radsmplr(ip, midpt, r(i), 1);
    smpl(isnan(smpl)) = 0;
    rp(i) = mean(smpl);
end

%% re-map
% generate the meshgrid for radial lookup
[vx, vy] = meshgrid(1:imSz(2), 1:imSz(1));
vx = vx - midpt(1);
vy = vy - midpt(2);
% calcualte the radius matrix
vr = hypot(vx, vy);
% mapping, out of range data are sampled as
I = interp1(r, rp, vr, 'spline', 0);

end
