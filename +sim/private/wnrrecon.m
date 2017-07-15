function J = wnrrecon(Fopt, imSz, popt, pr, parms)
%WNRRECON Combined images through a generalized Wiener filter.
%   
%   TBA

persistent A;

% create the apodization function if not exists
if isempty(A)
    A = filter.tukeywin2(imSz, parms.ApodizeRatio);
    A = single(A);
    
%     figure('Name', 'Apodization Function', 'NumberTitle', 'off');
%     imagesc(A);
%         axis image;
end

% extract frequent use parameters
nOri = parms.Orientations;
nPhase = parms.Phases;
w = parms.WienerConstant;

% interpolated size
rSz = parms.RetrievalInterpRatio*imSz;

%% OTF power spectrum
% FT
Iotf = psf2otf(parms.PSF, imSz);
Iotf = single(Iotf);
Iotf = ifftshift(Iotf);
% power spectrum
Potf = Iotf .* conj(Iotf);

%% estimate noise power spectrum
% cutoff frequency
f = 2*parms.NA / parms.Wavelength;
% chip size (real world scale) in nm
cSz = imSz .* parms.PixelSize;
% cutoff radius, short/long axis
r = f ./ (1 ./ cSz);

% generate the mask for region outside OTF cutoff frequency
[vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
vx = vx - floor(imSz(1)/2) + 1;
vy = vy - floor(imSz(2)/2) + 1;

m = ones(imSz, 'single');
m(((vx/r(1)).^2 + (vy/r(2)).^2) <= 1) = 0;
nElem = sum(m(:));

%% Wiener denominator
C = zeros([imSz, nPhase, nOri], 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % extract volume
        T = Fopt(:, :, iPhase, iOri);
        
        T = (T ./ Iotf);
        
        % calculate the average noise spectrum
        Pn = T .* m;
        Pn = Pn .* conj(Pn);
        Sn = sum(Pn(:)) / nElem;
        
        % mean power spectrum
        T = mean(T(:));
        T = (abs(T))^2;
        
        C(:, :, iPhase, iOri) = T * Potf / Sn;
    end
end
SC = sum(reshape(C, [imSz, nPhase*nOri]), 3);

% calculate the coefficients
C = C ./ (w+SC);

%% sum by Weiner filter
J = zeros(rSz, 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % extract volume
        T = Fopt(:, :, iPhase, iOri);
        
        %% multiply the filter function
        T = (T ./ Iotf) .* C(:, :, iPhase, iOri);
        
        % upsampling to perform FT interpolation in real space
        li = floor((rSz-imSz)/2)+1;
        ui = li+imSz-1;
        U = zeros(rSz, 'single');
        U(li(1):ui(1), li(2):ui(2), :) = T;

        %% apply relative phase shift
        if iPhase > 1
            U = fftshift(ifft2(ifftshift(U)));
            U = U .* pr(:, :, iPhase-1);
            U = fftshift(fft2(ifftshift(U)));
        end
        
        %% sum
        J = J + U;
    end
end

% convert back to the real space
J = fftshift(ifft2(ifftshift(J), 'symmetric'));

end
