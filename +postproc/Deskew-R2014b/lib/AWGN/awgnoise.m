function y = awgnoise(x, snr)
%AWIGNOISE Add averaged white gaussian noise to the signal.
%   
%   Y = AWGNOISE(X, SNR) adds AWG noise to signal X to generate a resulting
%   signal vector Y of specified SNR in dB.

% Set the random generator seed to default.
rng('default'); rng('shuffle');

L = length(x);

% Convert SNR to linear scale.
SNR = 10^(snr/10);
% Calculate the actual symbol energy.
Esym = sum(abs(x).^2)/(L);
% Find the noise spectral density.
N0 = Esym/SNR;

if isreal(x)
    % Standard deviation for AWGN when X is real.
    noiseSigma = sqrt(N0);
    n = noiseSigma*randn(1,L);
else
    % Standard deviation for AWGN when X is complex.
    noiseSigma = sqrt(N0/2);
    n = noiseSigma*(randn(1,L)+1i*randn(1,L));
end    

% Mixed signal. 
y = x + n;

end
