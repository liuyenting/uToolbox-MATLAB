function [res, varargout] = frcc2res(frcFrq, frcAvg)
%FRCC2RES Derive the resolution from an FRC curve.
%
%   RES = FRCC2RES(FRCFRQ, FRCAVG) finds the intersection between the 
%   averaged FRC curve and the 0.5-bit curve, deducted by information 
%   theory. RES is the estimated resolution in terms of the inverse unit of
%   spatial frequency used by FRCFRQ.
%
%   [RES, FRCTHR] = FRCC2RES(...) provides the additional threshold curve
%   through the output variable FRCTHR.

% the threshold curve
frcThr = (exp(1) / (6*pi)) * ones(size(frcFrq));

% intersection
lineThr = [frcFrq; frcThr];
lineAvg = [frcFrq; frcAvg];
P = interx(lineAvg, lineThr);

% use the first point only
P = P(1, :);

% lookup the spatial frequency
res = 1/P(1);

if nargout == 2
    varargout{1} = frcThr;
end

end
