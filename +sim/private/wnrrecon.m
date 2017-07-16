function J = wnrrecon(F, imSz, popt, prs, pr, parms)
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

%% create the initial phase shifts
popt = reshape(popt, [1, (nPhase-1)/2 * nOri]);
popt = exp(1i * [+popt; -popt]);
% reshape the matrix to match F layout
popt = popt(:);
popt = reshape(popt, [nPhase-1, nOri]);

%% calculate the filter functions
% system OTF
OTFsys = psf2otf(parms.PSF, imSz);
OTFsys = single(OTFsys);
% revert to natural order
OTFsys = ifftshift(OTFsys);

%% apply initial phase shift
Iotf = zeros([imSz, nPhase, nOri], 'single');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        T = OTFsys;

        % apply initial phase shift
        if iPhase > 1
            T = T .* popt(iPhase-1, iOri);
        end
        
        Iotf(:, :, iPhase, iOri) = T;
    end
end

%% generate Wiener coefficients
% nominators and (partial) denominators of filter functions
FFn = zeros([imSz, nPhase, nOri], 'single');

% padding offsets
li = floor((rSz-imSz)/2)+1;
ui = li+imSz-1;


% fill the nominator
h = figure('Name', 'Shifted Apodization Function', 'NumberTitle', 'off');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % shift the apodization function
        As = A;
        if iPhase > 1
            % apodization function buffer
            Ap = zeros(rSz, 'single');

            % pad
            Ap(li(1):ui(1), li(2):ui(2)) = As;
            
            Ap = fftshift(ifft2(ifftshift(Ap)));
            Ap = Ap .* conj(pr(:, :, iPhase-1, iOri));
            Ap = fftshift(fft2(ifftshift(Ap)));
            
            % crop
            As = Ap(li(1):ui(1), li(2):ui(2));
        end

        % create the nominator
        FFn(:, :, iPhase, iOri) = conj(Iotf(:, :, iPhase, iOri)) .* As;
        
        %% plotter
        figure(h);
        subplot(1, nPhase, iPhase);
            imagesc(abs(As).^0.1);
            axis image;
        if iPhase == 1
            title('m_0');
        else 
            if mod(iPhase, 2) == 0
                s = '-';
            else
                s = '+';
            end
            t = sprintf('m_%d^%c', floor(iPhase/2), s);
            title(t);
        end
    end
end

% buffer space for upsampled frequency components
Fp = zeros([rSz, nPhase, nOri], 'single');

h = figure('Name', 'Filtered Components', 'NumberTitle', 'off');
for iOri = 1:nOri
    for iPhase = 1:nPhase
        % calculate the denominator
        FFd = zeros(imSz, 'single');
        for iiPhase = 2:nPhase
            % no need to calculate itself (subtraction yields 0)
            if iPhase == iiPhase 
                continue;
            end
            
            % phase shift differences (m0 doesn't have original offset)
            D = prs(iiPhase-1, iOri);
            if iPhase > 1
                D = D .* conj(prs(iPhase-1, iOri));
            end
            FFd = FFd + abs(Iotf(:, :, iPhase, iOri) .* D).^2;
        end
        % include the Wiener constant
        FFd = FFd + w;
        
        % save the filtered element (and zero-padded)
        Fp(li(1):ui(1), li(2):ui(2), iPhase, iOri) = ...
            (FFn(:, :, iPhase, iOri) ./ FFd) .* F(:, :, iPhase, iOri);
        
        %% plotter
        figure(h);
        subplot(1, nPhase, iPhase);
            imagesc(abs(Fp(li(1):ui(1), li(2):ui(2), iPhase, iOri)).^0.1);
            axis image;
        if iPhase == 1
            title('m_0');
        else 
            if mod(iPhase, 2) == 0
                s = '-';
            else
                s = '+';
            end
            t = sprintf('m_%d^%c', floor(iPhase/2), s);
            title(t);
        end
    end
end

%% apply the relative phase shifts
T = fftshift(ifft2(ifftshift(Fp(:, :, 2:end, :))));
T = T .* conj(pr);
Fp(:, :, 2:end, :) = fftshift(fft2(ifftshift(T)));

%% plotter
h = figure('Name', 'Relative Phase Shifts (applied)', 'NumberTitle', 'off');
for iPhase = 1:nPhase
    figure(h);
    subplot(1, nPhase, iPhase);
        imagesc(abs(Fp(:, :, iPhase, :)).^0.1);
        axis image;
    if iPhase == 1
        title('m_0');
    else 
        if mod(iPhase, 2) == 0
            s = '-';
        else
            s = '+';
        end
        t = sprintf('m_%d^%c', floor(iPhase/2), s);
        title(t);
    end
end

%% add together and transform back to yield the answer
Fp = reshape(Fp, [rSz, nPhase*nOri]);
J = sum(Fp, 3);

figure('Name', 'Merged', 'NumberTitle', 'off');
subplot(1, 2, 1);
imagesc(abs(J).^0.1);
    axis image;
    title('K-Space');

J = fftshift(ifft2(ifftshift(J), 'symmetric'));

subplot(1, 2, 2);
imagesc(J);
    axis image;
    title('R-Space');

disp('=== HELLO WORLD ===');

% %% OTF power spectrum
% % FT
% Iotf = psf2otf(parms.PSF, imSz);
% Iotf = single(Iotf);
% Iotf = ifftshift(Iotf);
% % power spectrum
% Potf = Iotf .* conj(Iotf);
% 
% %% estimate noise power spectrum
% % cutoff frequency
% f = 2*parms.NA / parms.Wavelength;
% % chip size (real world scale) in nm
% cSz = imSz .* parms.PixelSize;
% % cutoff radius, short/long axis
% r = f ./ (1 ./ cSz);
% 
% % generate the mask for region outside OTF cutoff frequency
% [vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
% vx = vx - floor(imSz(1)/2) + 1;
% vy = vy - floor(imSz(2)/2) + 1;
% 
% m = ones(imSz, 'single');
% m(((vx/r(1)).^2 + (vy/r(2)).^2) <= 1) = 0;
% nElem = sum(m(:));
% 
% %% Wiener denominator
% C = zeros([imSz, nPhase, nOri], 'single');
% for iOri = 1:nOri
%     for iPhase = 1:nPhase
%         % extract volume
%         T = Fopt(:, :, iPhase, iOri);
%         
%         T = (T ./ Iotf);
%         
%         % calculate the average noise spectrum
%         Pn = T .* m;
%         Pn = Pn .* conj(Pn);
%         Sn = sum(Pn(:)) / nElem;
%         
%         % mean power spectrum
%         T = mean(T(:));
%         T = (abs(T))^2;
%         
%         C(:, :, iPhase, iOri) = T * Potf / Sn;
%     end
% end
% SC = sum(reshape(C, [imSz, nPhase*nOri]), 3);
% 
% % calculate the coefficients
% C = C ./ (w+SC);
% 
% %% sum by Weiner filter
% J = zeros(rSz, 'single');
% for iOri = 1:nOri
%     for iPhase = 1:nPhase
%         % extract volume
%         T = Fopt(:, :, iPhase, iOri);
%         
%         %% multiply the filter function
%         T = (T ./ Iotf) .* C(:, :, iPhase, iOri);
%         
%         % upsampling to perform FT interpolation in real space
%         li = floor((rSz-imSz)/2)+1;
%         ui = li+imSz-1;
%         U = zeros(rSz, 'single');
%         U(li(1):ui(1), li(2):ui(2), :) = T;
% 
%         %% apply relative phase shift
%         if iPhase > 1
%             U = fftshift(ifft2(ifftshift(U)));
%             U = U .* pr(:, :, iPhase-1);
%             U = fftshift(fft2(ifftshift(U)));
%         end
%         
%         %% sum
%         J = J + U;
%     end
% end

% %% simple sum
% popt = reshape(popt, [1, (nPhase-1)/2 * nOri]);
% popt = exp(1i * [+popt; -popt]);
% % flatten the array for linear duplication
% popt = popt(:);
% % apply initial phase shift
% F(:, :, 2:end, :) = F(:, :, 2:end, :) .* reshape(popt, [1, 1, nPhase-1, nOri]);
% 
% Fp = zeros([rSz, nPhase, nOri], 'single');
% % upsampling to perform FT interpolation in real space
% li = floor((rSz-imSz)/2)+1;
% ui = li+imSz-1;
% Fp(li(1):ui(1), li(2):ui(2), :, :) = F;
% 
% % back to time domain
% Fp(:, :, 2:end, :) = fftshift(ifft2(ifftshift(Fp(:, :, 2:end, :))));
% 
% % add relative phase shift deduced from kp values (imaginary number in the
% % time domain)
% Fp(:, :, 2:end, :) = Fp(:, :, 2:end, :) .* pr;
% 
% Fp(:, :, 2:end, :) = fftshift(fft2(ifftshift(Fp(:, :, 2:end, :))));
% 
% % sum all the orientations and phases
% Fp = reshape(Fp, [rSz, nPhase*nOri]);
% J = sum(Fp, 3);
% 
% % convert back to the real space
% J = fftshift(ifft2(ifftshift(J), 'symmetric'));

end
