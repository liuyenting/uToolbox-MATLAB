function M = spectramat(np, I0, I1)
%SPECTRAMAT Generate spectra matrix for the linear system.
%   
%   M = SPECTRAMAT(NP, I0, I1) generates the coefficient matrix using
%   number of phases and beam intensities. Use alongside with the phase
%   solver.
%
%   See also: TBA

%TODO modify to adpat to N-D SIM

% create the phase array
P = linspace(0, 2*pi, np+1);
P = P(1:end-1).';

% 1st order uses dummy exp(0) to avoid horizontal concatenation
%   m0, m1-, m1+, m2-, m2+
M = [ ...
    (I0+2*I1) * exp(0 * P), ...
    sqrt(I0*I1) * exp(1i * P), ...
    sqrt(I0*I1) * exp(-1i * P), ...
    I1 * exp(1i * 2*P), ...
    I1 * exp(-1i * 2*P) ...
];

% divide a constant factor which scale the intensity of the captured image
% (assusmed to be) 1
M = M / (I0+2*I1);

% type cast to single precision in order to avoid unnecessary type cast in
% rest of the pipeline
M = single(M);

end

