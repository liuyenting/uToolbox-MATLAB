function y = awgnoise(x, snr)
%AWIGNOISE Add averaged white gaussian noise to the signal.
%   
%   Y = AWGNOISE(X, SNR) adds AWG noise to signal X to generate a resulting
%   signal vector Y of specified SNR in dB.
%
%   Reference
%   ---------
%   How to generate AWGN noise in Matlab/Octave (without using in-built
%   awgn function). Mathuranathan Viswanathan. June 16, 2015.

%% seeding based on current time
rng('shuffle');

%% pre-calculate
L = length(x);

% convert SNR to linear scale
snr = 10^(snr/10);
% calculate the actual symbol energy
Esym = sum(abs(x).^2)/(L);
% find the noise spectral density
n0 = Esym/snr;

%% generate
% SD for the noise (different case in real/complex data)
if isreal(x)
    sigma = sqrt(n0);
    n = sigma*randn(1,L);
else
    sigma = sqrt(n0/2);
    n = sigma*(randn(1,L)+1i*randn(1,L));
end    

%% apply the noise
y = x + n;

end
