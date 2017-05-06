function [res, frcThr] = frcc2res(frcFrq, frcAvg)
%FRCC2RES Derive the resolution from an FRC curve.
%
%   TBA

% the threshold curve
frcThr = (exp(1) / (6*pi)) * ones(size(frcFrq));

% intersection
lineThr = [frcFrq; frcThr];
lineAvg = [frcFrq; frcAvg];
ind = interx(lineThr, lineAvg);

res = 1/ind(1);

end
