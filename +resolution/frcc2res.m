function [res, frcThr] = frcc2res(frcFrq, frcAvg)
%FRCC2RES Derive the resolution from an FRC curve.
%
%   TBA

% the threshold curve
frcThr = (exp(1) / (6*pi)) * ones(size(frcFrq));

% intersection
lineThr = [frcFrq; frcThr];
lineAvg = [frcFrq; frcAvg];
P = interx(lineThr, lineAvg);

% use the first point only
P = P(1, :);

% lookup the spatial frequency
res = 1/P(1);

end
