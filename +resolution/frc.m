function varargout = frc(A, B, tsz)
%FRC Calculate Fourier ring correlation values.
%
%   COR = FRC(A, B)
%   COR = FRC(A, B, TSZ) calculates the correlation values COR along the
%   axis of the incoming images A and B. TSZ indicates fractions of the
%   image that will be damped during the masking. TSZ is default to 8,
%   meaning that 1/8 * 2 = 2/8 of the images, along the border line, will
%   get damped. 

npx = size(A);
if size(B) ~= npx
    error('resolution:frc', ...
          'Image sizes do not match.');
end

if nargin == 2
    tsz = 8;
end

% generate Tukey window
mask = tukeywin2(npx, tsz);

% mask the binned images
A = A .* mask;
B = B .* mask;

% calculate the FT result
FA = fftshift(fft2(A));
FB = fftshift(fft2(B));

% numerator
num = radialsum(FA.*conj(FB));
num = real(num);

% denominator
SA = radialsum(abs(FA).^2);
SB = radialsum(abs(FB).^2);
den = sqrt(abs(SA.*SB));

% divided result
cor = num ./ den;
% remove NaN
cor(isnan(cor)) = 0;

% set the output
varargout{1} = cor;
if nargout == 2
    varargout{2} = num;
end

% DEBUG
%printft(A, B, FA, FB);

end

function printft(A, B, FA, FB) %#ok<DEFNU>
%PRINTSMPLTRACE Print the sampled trace on source array.

figure('Name', '[DEBUG] FFT Result', 'NumberTitle', 'off'); 

subplot(2, 2, 1);
imagesc(A);
axis image;
title('A');

subplot(2, 2, 2);
imagesc(B);
axis image;
title('B');

FA = 100*log(1+abs(fftshift(FA)));
subplot(2, 2, 3);
imagesc(FA);
axis image;
title('F(A)');

FB = 100*log(1+abs(fftshift(FB)));
subplot(2, 2, 4);
imagesc(FB);
axis image;
title('F(B)');

end
