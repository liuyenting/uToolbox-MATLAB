function M = spectramat(np, I0, I1)
%SPECTRAMAT Generate spectra matrix for the linear system.
%   
%   M = SPECTRAMAT(NP, I0, I1) generates the coefficient matrix using
%   number of phases and beam intensities. Use alongside with the phase
%   solver.
%
%   See also: TBA

% create the phase array
P = linspace(0, 2*pi, np+1);
P = P(1:end-1).';

% 1st order uses dummy exp(0) to avoid horizontal concatenation
%   m0, m-1, m+1, m-2, m+2
M = [ ...
    (I0+2*I1) * exp(0 * P), ...
    sqrt(I0*I1) * exp(-1i * P), ...
    sqrt(I0*I1) * exp(1i * P), ...
    I1 * exp(-1i * 2*P), ...
    I1 * exp(1i * 2*P) ...
];

end

