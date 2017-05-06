function [res, frc_thr] = frcc2res(frc_frq, frc_avg)
%FRCC2RES Derive the resolution from an FRC curve.
%
%   TBA

% the threshold curve
frc_thr = (exp(1) / (6*pi)) * ones(size(frc_frq));

% intersection
dat_thr = [frc_frq; frc_thr];
dat_avg = [frc_frq; frc_avg];
ind = interx(dat_thr.', dat_avg.');

res = 1/ind(1);

end
