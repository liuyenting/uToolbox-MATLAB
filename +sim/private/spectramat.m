function M = spectramat(np, varargin)
%SPECTRAMAT Generate spectra matrix for the linear system.
%   
%   M = SPECTRAMAT(NP, I0, I1) generates the coefficient matrix using
%   number of phases and beam intensities. Use alongside with the phase
%   solver.
%
%   See also: TBA

%% parameters
% identify the type of SIM
%   * 2-D SIM requires 3 phases
%   * 3-D SIM requires 5 phases
if np == 3
    nd = 2;
elseif np == 5
    nd = 3;
else
    error(generatemsgid('InvalidPhases'), ...
          'Invalid amount of phases provided.');
end

% set the intensity component
if nd == 3
    if nargin > 1
        I0 = varargin{1};
    else
        I0 = 1;
    end
    if nargin > 2
        I1 = varargin{2};
    end
end

if nargin > 3 || ((nd == 2) && (nargin > 2))
    warning(generatemsgid('ExcessVarIn'), ...
            'Excessive input parameters are ignored.');
end

%% process
% create the phase array
P = linspace(0, 2*pi, np+1);
P = P(1:end-1).';

% 1st order uses dummy exp(0) to avoid horizontal concatenation
%   m0, m1-, m1+, m2-, m2+
% 
% matrix M is divided by a constant factor which scale (and assume) the 
% intensity of the captured image to 1
if nd == 2
    M = [exp(0 * P), exp(1i * P), exp(-1i * P)];
elseif nd == 3 
    M = [ ...
        (I0+2*I1) * exp(0 * P), ...
        sqrt(I0*I1) * exp(1i * P), ...
        sqrt(I0*I1) * exp(-1i * P), ...
        I1 * exp(1i * 2*P), ...
        I1 * exp(-1i * 2*P) ...
    ];
    % divide a constant factor which scale the intensity of the captured 
    % image (assusmed to be) 1
    M = M / (I0+2*I1);
end

% type cast to single precision in order to avoid unnecessary type cast in
% rest of the pipeline
M = single(M);

end

